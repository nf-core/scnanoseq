process COMBINE_SEURAT_STATS {
    label 'process_low'

    //TODO: figure out env.

    input:
    path seurat_stats

    output:
    path "*.tsv"        , emit: combined_stats
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    combine_files.sh \\
    $args \\
    -i ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_seurat_stats: v1.0
    END_VERSIONS
    """
}
