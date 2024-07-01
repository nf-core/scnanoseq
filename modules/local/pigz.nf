process PIGZ {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4':
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(unzipped_file)
    val addl_prefix

    output:
    tuple val(meta), path("*.gz"), emit: archive
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${addl_prefix}"
    """
    cat ${unzipped_file} > ${prefix}.fastq
    pigz \\
        $args \\
        -f \\
        -p $task.cpus \\
        ${prefix}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz -V 2>&1) | sed 's/pigz //g' )
    END_VERSIONS
    """
}
