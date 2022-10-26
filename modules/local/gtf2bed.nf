process GTF2BED {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '\$1 !~ /^#/ { print \$1, \$4-1, \$5 }' > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort_gtf: 1.0
    END_VERSIONS
    """
}
