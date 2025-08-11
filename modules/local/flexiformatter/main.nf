process FLEXIFORMATTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiformatter%3A1.0.3--pyhdfd78af_0':
        'https://quay.io/biocontainers/flexiformatter:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    val bam_index_extension

    output:
    tuple val(meta), path("*_tagged.bam")       , emit: bam
    tuple val(meta), path("*_tagged.bam.bai")   , emit: bai
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_tagged"
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.bam"
    def bam_output = " | samtools sort -@ ${task.cpus-1} -o ${bam_index} ${args2}"

    """
    flexiformatter \\
        ${bam} \\
        ${args} \\
        ${bam_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiformatter: \$(echo \$(flexiformatter --version) |& sed 's/flexi_formatter version: //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiformatter: \$(flexiformatter --version |& sed 's/flexi_formatter version: //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
