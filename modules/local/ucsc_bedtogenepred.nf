process UCSC_BEDTOGENEPRED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ucsc-bedtogenepred=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtogenepred:377--h0b8a92a_2' :
        'biocontainers/ucsc-bedtogenepred:377--h0b8a92a_2' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.pred"), emit: pred
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    awk '{print \$0"\t."}' ${bed} > ${bed}.processed

    bedToGenePred \\
        ${bed}.processed \\
        ${prefix}.pred
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
