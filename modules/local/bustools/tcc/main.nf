process BUSTOOLS_TCC {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bustools=0.45.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bustools:0.45.1--h6f0a7f7_0' :
        'biocontainers/bustools:0.45.1--h6f0a7f7_0' }"

    input:
    tuple val(meta), path(bus), path(ecmap), path(txnames), path(known_barcodes)
    tuple val(meta2), path(t2g)
    val correct

    output:
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.mtx")         , emit: tcc_mtx
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.ec.txt")      , emit: tcc_ec
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.barcodes.txt"), emit: barcodes
    tuple val(meta), path("gene/matrix.mtx.gz")                        , emit: gene_mtx
    tuple val(meta), path("gene/features.tsv.gz")                      , emit: gene_features
    tuple val(meta), path("gene/barcodes.tsv.gz")                      , emit: gene_barcodes
    path "versions.yml"                                                , emit: versions_bustools_tcc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args reaches the equivalence-class count only and ext.args2 the gene
    // count. They are kept apart because the two counts must differ on
    // --multimapping: the TCC count needs it, the gene count must not have it,
    // or gene-ambiguous UMIs come back as fractions in a matrix declared integer.
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def memory = task.memory.toGiga() - 1
    // Correcting against the known-barcode list is also what discards the reads
    // that matched no known barcode, which is what makes the counts called-cells
    // only. The all-droplet pass skips it: its barcodes come from XB, which is
    // already the final droplet identity, and every one of them has to survive.
    def correct_step = !correct ? "mv sorted.bus counted.bus" : """
    # Reduce the flexiplex known-barcode list to one barcode per line. The
    # barcodes are already corrected against the whitelist, so this on-list
    # stops bustools inventing one of its own.
    awk '{print \$1}' ${known_barcodes} \\
        | grep -E '^[ACGTN]+\$' \\
        | sort -u \\
        > onlist.txt

    bustools correct \\
        -o corrected.bus \\
        -w onlist.txt \\
        sorted.bus

    bustools sort \\
        -o counted.bus \\
        -T tmp \\
        -t ${task.cpus} \\
        -m ${memory}G \\
        corrected.bus
    """
    """
    mkdir -p tmp counts_unfiltered genecounts

    bustools sort \\
        -o sorted.bus \\
        -T tmp \\
        -t ${task.cpus} \\
        -m ${memory}G \\
        ${bus}

    ${correct_step}

    # The sorted (and, for called cells, barcode-corrected) BUS file is counted
    # twice. Neither count sees a deduplicated file: each one collapses UMIs per
    # gene itself with --umi-gene, and because only the first keeps
    # gene-ambiguous records the two do not collapse the same set.
    #
    # First at the equivalence class level, which is what kallisto quant-tcc
    # consumes for the transcript EM: no --genecounts, and --multimapping so
    # that records compatible with more than one gene reach the EM rather than
    # being dropped.
    bustools count \\
        -o counts_unfiltered/cells_x_tcc \\
        -g ${t2g} \\
        -e ${ecmap} \\
        -t ${txnames} \\
        --umi-gene \\
        --multimapping \\
        ${args} \\
        counted.bus

    # Then at the gene level for the published gene matrix. Without
    # --multimapping a UMI compatible with more than one gene is discarded
    # rather than split, so every entry is a whole UMI. This is the standard
    # kb-python gene matrix and the one to feed to CellBender or any other tool
    # with a count likelihood. The quant-tcc EM gene abundances, which
    # distribute those UMIs instead, are fractional and stay available in
    # quant/matrix.abundance.gene.mtx.
    bustools count \\
        -o genecounts/cells_x_genes \\
        -g ${t2g} \\
        -e ${ecmap} \\
        -t ${txnames} \\
        --genecounts \\
        --umi-gene \\
        ${args2} \\
        counted.bus

    # bustools names genes by gene_id alone. CellBender's MEX reader indexes
    # three columns whenever features.tsv.gz exists, and Read10X treats such a
    # directory as Cell Ranger v3, so add gene_name from the t2g (column 3)
    # and the 10x feature type. Genes without a name fall back to their id.
    awk -F'\\t' -v OFS='\\t' '
        NR == FNR { name[\$2] = \$3; next }
        { print \$1, (\$1 in name ? name[\$1] : \$1), "Gene Expression" }
    ' ${t2g} genecounts/cells_x_genes.genes.txt > genecounts/features.tsv

    # bustools writes cells x genes; Read10X and CellBender want features x
    # cells. --integer checks every value is a whole UMI before declaring the
    # matrix integer, so a stray fraction fails the task instead of being
    # silently truncated downstream.
    mtx_transpose_to_mex.sh \\
        --mtx genecounts/cells_x_genes.mtx \\
        --features genecounts/features.tsv \\
        --barcodes genecounts/cells_x_genes.barcodes.txt \\
        --outdir gene \\
        --integer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version | sed 's/^bustools, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p counts_unfiltered gene
    touch counts_unfiltered/cells_x_tcc.mtx
    touch counts_unfiltered/cells_x_tcc.ec.txt
    touch counts_unfiltered/cells_x_tcc.barcodes.txt
    echo "" | gzip -c > gene/matrix.mtx.gz
    echo "" | gzip -c > gene/features.tsv.gz
    echo "" | gzip -c > gene/barcodes.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version | sed 's/^bustools, version //')
    END_VERSIONS
    """
}
