process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::minimap2=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.24--h5bf99c6_0':
        'biocontainers/minimap2:2.24--h5bf99c6_0' }"

    input:
    tuple val(meta), path(fastq)
    path bed
    path reference

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def junctions = "--junc-bed ${bed}"
    """
    minimap2 \\
        $args \\
        $junctions \\
        -t $task.cpus \\
        $reference \\
        $fastq > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
