process GET_GTF_FEATURES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gtf)
    val feature

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk 'OFS="\\t", \$1 !~ /^#/ {if (\$3 == "exon") print \$0}' ${gtf} > ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_gtf_features: 1.0
    END_VERSIONS
    """
}
