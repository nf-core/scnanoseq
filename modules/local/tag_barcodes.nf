process TAG_BARCODES {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pysam=0.19.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.19.1--py310hff46b53_1' :
        'biocontainers/pysam:0.19.1--py310hff46b53_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(corrected_bc_info)

    output:
    tuple val(meta), path("*.tagged.bam"), emit: tagged_bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tag_barcodes.py \\
        -b ${bam} \\
        -i ${corrected_bc_info} \\
        -o ${prefix}.tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
