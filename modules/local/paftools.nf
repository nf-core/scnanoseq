process PAFTOOLS {
    label 'process_low'
    // TODO: decide what's best tag after more gtf inputs have been added?

    conda (params.enable_conda ? "bioconda::minimap2=2.24" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.24--h5bf99c6_0':
        'quay.io/biocontainers/minimap2:2.24--h5bf99c6_0' }"

    // TODO: *** once intron method 1/2 gets added, add conditionals to input below ***
    input:
    path gtf

    output:
    path "*.bed", emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    paftools.js gff2bed \\
        $args \\
        $gtf > ${gtf.baseName}.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

