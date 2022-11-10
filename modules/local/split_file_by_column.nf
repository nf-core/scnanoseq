process SPLIT_FILE_BY_COLUMN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(unsplit_file)
    val split_amount

    output:
    tuple val(meta), path("*.tsv"), emit: split_files
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_file_by_column.sh -i ${unsplit_file} -s $split_amount

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split_file_by_col: v1.0 
    END_VERSIONS
    """
}
