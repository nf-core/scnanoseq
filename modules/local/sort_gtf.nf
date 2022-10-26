process SORT_GTF {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.sorted.gtf"), emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '\$1 !~ /^#/ { print \$0 }' | sort -k1,1n -k4,4n -k5,5n > ${prefix}.sorted.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort_gtf: 1.0
    END_VERSIONS

    """
}
