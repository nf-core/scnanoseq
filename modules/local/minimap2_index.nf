process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::minimap2=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.24--h5bf99c6_0':
        'biocontainers/minimap2:2.24--h5bf99c6_0' }"

    input:
    path fasta
    path bed

    output:
    path "*.mmi"        , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // TODO: see if there is a better way of including additional
    // input (e.g.: bed / junctions), so we can use the module in nf-core rather than local

    script:
    def args = task.ext.args ?: ''
    def junctions = "--junc-bed ${bed}"
    """
    minimap2 \\
        $args \\
        $junctions \\
        -t $task.cpus \\
        -d ${fasta}.mmi \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
