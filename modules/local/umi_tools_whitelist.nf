
process UMI_TOOLS_WHITELIST {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38hbff2b2d_1':
        'quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: whitelist
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    umi_tools \\
        whitelist \\
        --log2stderr \\
        --stdin=${reads} \\
        ${args} > ${whitelist} 2> /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitoolswhitelist: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
