process PROWLERTRIMMER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "jdoe062894::prowlertrimmer" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    gunzip -c $reads > ${prefix}.fastq
    python3 \$(which TrimmerLarge.py) $args -f ${prefix}.fastq
    gzip ${prefix}TrimLT-U0-S${params.min_q_score}W100L${params.min_length}R0.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prowlertrimmer: \$(echo "1.0")
    END_VERSIONS
    """
}
