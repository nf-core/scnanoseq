process CORRECT_COUNTS_MATRIX {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pandas=1.5.1" : null)
    container "docker.io/biocontainers/pandas:1.5.1_cv1" // from PR: https://github.com/BioContainers/containers/pull/504

    input:
    tuple val(meta), path(counts_matrix)

    output:
    tuple val(meta), path("*.corrected.tsv"), emit: corrected_counts_matrix
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    FILE_PREFIX=${prefix}.\$(basename ${counts_matrix} | cut -f2 -d'.')


    correct_counts_matrix.py \\
        -i $counts_matrix \\
        -o \${FILE_PREFIX}.corrected.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
