process PROWLERTRIMMER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "jdoe062894::prowlertrimmer" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    gunzip -c $reads > ${prefix}.fastq
    python3 \$(which TrimmerLarge.py) $args -f ${prefix}.fastq -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prowlertrimmer: \$(echo "1.0")
    END_VERSIONS
    """
}
