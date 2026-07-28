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

    output:
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.mtx")         , emit: tcc_mtx
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.ec.txt")      , emit: tcc_ec
    tuple val(meta), path("counts_unfiltered/cells_x_tcc.barcodes.txt"), emit: barcodes
    path "versions.yml"                                                , emit: versions_bustools_tcc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory.toGiga() - 1
    """
    # Reduce the flexiplex known-barcode list to one barcode per line. The
    # barcodes are already corrected against the whitelist, so this on-list
    # stops bustools inventing one of its own.
    awk '{print \$1}' ${known_barcodes} \\
        | grep -E '^[ACGTN]+\$' \\
        | sort -u \\
        > onlist.txt

    mkdir -p tmp counts_unfiltered

    bustools sort \\
        -o sorted.bus \\
        -T tmp \\
        -t ${task.cpus} \\
        -m ${memory}G \\
        ${bus}

    bustools correct \\
        -o corrected.bus \\
        -w onlist.txt \\
        sorted.bus

    bustools sort \\
        -o corrected.sorted.bus \\
        -T tmp \\
        -t ${task.cpus} \\
        -m ${memory}G \\
        corrected.bus

    # Omitting --genecounts keeps this at the equivalence class level, which is
    # what kallisto quant-tcc consumes. --umi-gene deduplicates umis.
    bustools count \\
        -o counts_unfiltered/cells_x_tcc \\
        -g ${t2g} \\
        -e ${ecmap} \\
        -t ${txnames} \\
        --umi-gene \\
        --multimapping \\
        ${args} \\
        corrected.sorted.bus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version | sed 's/^bustools, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p counts_unfiltered
    touch counts_unfiltered/cells_x_tcc.mtx
    touch counts_unfiltered/cells_x_tcc.ec.txt
    touch counts_unfiltered/cells_x_tcc.barcodes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version | sed 's/^bustools, version //')
    END_VERSIONS
    """
}
