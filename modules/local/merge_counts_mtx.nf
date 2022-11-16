process MERGE_COUNTS_MTX {
    tag "$meta.id"
    label 'process_low'

    conda ("conda-forge::r-tidyverse=1.3.1 conda-forge::r-optparse")

    input:
    tuple val(meta), path(mtx_1), path(mtx_2)

    output:
    tuple val(meta), path("*.merged.tsv"), emit: merged_mtx
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_counts.R \\
        $args \\
        -i ${mtx_1} \\
        -j ${mtx_2} \\
        -o ${prefix}.merged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_counts_mtx: 1.0
    END_VERSIONS
    """
}
