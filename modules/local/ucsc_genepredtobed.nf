process UCSC_GENEPREDTOBED {
    tag "$genepred"
    label 'process_low'

    conda "bioconda::ucsc-genepredtobed=447"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-genepredtobed:447--h954228d_0':
        'biocontainers/ucsc-genepredtobed:447--h954228d_0' }"

    input:
    path genepred

    output:
    path "*.bed"        , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    genePredToBed \\
        $args \\
        $genepred \\
        ${genepred.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
