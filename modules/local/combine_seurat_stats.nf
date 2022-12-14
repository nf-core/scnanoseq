process COMBINE_SEURAT_STATS {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

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
        combine_seurat_stats: 1.0
    END_VERSIONS
    """
}
