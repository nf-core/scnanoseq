process MERGE_FILE_BY_COLUMN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pandas=1.5.1" : null)
    container "docker.io/biocontainers/pandas:1.5.1_cv1" // from PR: https://github.com/BioContainers/containers/pull/504

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

    merge_files_by_column.py -o ${prefix}.merged.tsv -i \$(pwd) -s tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
