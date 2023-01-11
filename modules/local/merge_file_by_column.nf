process MERGE_FILE_BY_COLUMN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

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
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
        cut: \$(echo \$(cut --version) | sed 's/^.*cut (GNU coreutils) //; s/ .*//')
        paste: \$(echo \$(paste --version) | sed 's/^.*paste (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
