process UCSC_GENEPREDTOGTF {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-genepredtogtf:377--h0b8a92a_2' :
        'quay.io/biocontainers/ucsc-genepredtogtf:377--h0b8a92a_2' }"

    input:
    tuple val(meta), path(pred)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    genePredToGtf \\
        $pred \\
        ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
