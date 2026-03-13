process COMBINE_SEURAT_STATS {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

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
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch combined_seurat.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
