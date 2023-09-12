process UCSC_GENEPREDTOGTF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ucsc-genepredtogtf=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-genepredtogtf:377--h0b8a92a_4' :
        'biocontainers/ucsc-genepredtogtf:377--h0b8a92a_4' }"

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
    genePredToGtf file \\
        $pred \\
        ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
