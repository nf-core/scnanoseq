process KALLISTO_QUANTTCC {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kallisto=0.52.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.52.0--h13ff97a_0' :
        'biocontainers/kallisto:0.52.0--h13ff97a_0' }"

    input:
    tuple val(meta), path(tcc_mtx), path(tcc_ec), path(barcodes)
    tuple val(meta2), path(index)
    tuple val(meta3), path(t2g)

    output:
    tuple val(meta), path("transcript/matrix.mtx.gz")  , emit: transcript_mtx
    tuple val(meta), path("transcript/features.tsv.gz"), emit: transcript_features
    tuple val(meta), path("transcript/barcodes.tsv.gz"), emit: transcript_barcodes
    tuple val(meta), path("quant/*")                   , emit: quant
    path "versions.yml"                                , emit: versions_kallisto_quanttcc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto quant-tcc \\
        --long \\
        -o quant \\
        -i ${index} \\
        -e ${tcc_ec} \\
        -g ${t2g} \\
        -t ${task.cpus} \\
        ${args} \\
        ${tcc_mtx}

    mkdir -p transcript

    # quant-tcc writes cells x features; Read10X wants features x cells, so
    # transpose the header dimensions and every entry. Only the transcript
    # matrix is published in MEX form: the gene matrix comes from bustools
    # count --genecounts in BUSTOOLS_TCC, because the EM gene abundances here
    # are fractional. They remain available as quant/matrix.abundance.gene.mtx.
    grep '^%' quant/matrix.abundance.mtx > transcript/matrix.mtx
    grep -v '^%' quant/matrix.abundance.mtx | awk '{print \$2" "\$1" "\$3}' >> transcript/matrix.mtx

    cp quant/transcripts.txt transcript/features.tsv
    cp ${barcodes} transcript/barcodes.tsv

    gzip transcript/matrix.mtx transcript/features.tsv transcript/barcodes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p transcript quant
    echo "" | gzip -c > transcript/matrix.mtx.gz
    echo "" | gzip -c > transcript/features.tsv.gz
    echo "" | gzip -c > transcript/barcodes.tsv.gz
    touch quant/matrix.abundance.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """
}
