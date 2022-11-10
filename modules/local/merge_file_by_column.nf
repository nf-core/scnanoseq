process MERGE_FILE_BY_COLUMN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(split_files)

    output:
    tuple val(meta), path("*.merged.tsv"), emit: col_merged_file
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    merge_files_by_column.sh ${prefix}.merged.tsv $split_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_file_by_col: v1.0 
    END_VERSIONS
    """
}
