process GTF_TO_T2G {
    tag "gtf_to_t2g"
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.t2g.tsv"), emit: t2g
    path "versions.yml"               , emit: versions_gtf_to_t2g, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    gtf_to_t2g.py \\
        -i ${gtf} \\
        -o ${prefix}.t2g.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.t2g.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
