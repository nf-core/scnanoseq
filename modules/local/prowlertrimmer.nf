process PROWLERTRIMMER {
    tag "$meta.id"
    label 'process_low'

    conda "jdoe062894::prowlertrimmer"

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
    FILE_PREFIX=${prefix}
    if [ ${params.split_amount} -gt 0 ]; then
        IDX=\$(basename ${reads} | cut -f2 -d'.')
        FILE_PREFIX=\${FILE_PREFIX}.\${IDX}
    fi

    python3 \$(which TrimmerLarge.py) $args -f \${FILE_PREFIX}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prowlertrimmer: \$(echo "1.0")
    END_VERSIONS
    """
}
