process REFORMAT_WHITELIST {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(whitelist)

    output:
    tuple val(meta), path("*.bc_count.txt")        , emit: bc_list_counts
    tuple val(meta), path("*.bc_count.txt.bc_only"), emit: bc_list
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    reformat_whitelist.py \\
        -i ${whitelist} \\
        -o ${prefix}.bc_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
