process UCSC_GTFTOGENEPRED {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::ucsc-gtftogenepred=447"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:447--h954228d_0':
        'biocontainers/ucsc-gtftogenepred:447--h954228d_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    path "*.genepred"    , emit: genepred
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${gtf.baseName}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '447'
    """
    gtfToGenePred \\
        $args \\
        $gtf  \\
        ${prefix}.genepred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${gtf.baseName}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '447'
    """
    touch ${prefix}.genepred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
