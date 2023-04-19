process SPLIT_FILE_BY_COLUMN {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

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
        cut: \$(echo \$(cut --version) | sed 's/^.*cut (GNU coreutils) //; s/ .*//')
        paste: \$(echo \$(paste --version) | sed 's/^.*paste (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
