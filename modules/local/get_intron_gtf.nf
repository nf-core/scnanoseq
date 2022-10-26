process GET_INTRON_GTF {
    tag "$meta.id"
    label 'process_low'

    input:
    path chr_lengths
    path intergenic_regions
    path exon_regions

    output:
    tuple path("introns.gtf")     , emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    :q

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_intron_gtf: 1.0
    END_VERSIONS
    """
}
