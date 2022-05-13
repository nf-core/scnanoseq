process CREATE_REGEX {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'quay.io/biocontainers/python:3.8.3' }"

    output:
    path "regex_patterns.txt", emit: regex_pattern
    path "versions.yml"      , emit: versions

    when:
    params.identifier_pattern != null &&
    params.cell_barcode_lengths != null &&
    params.umi_lengths != null &&
    params.fixed_seqs != null

    script:
    def args = task.ext.args ?: ''
    
    """
    python create_regex.py -i ${params.identifier_pattern} \\
                           -c ${params.cell_barcode_lengths} \\
                           -u ${params.umi_lengths} \\
                           -f ${params.fixed_seqs} \\
                           -o regex_patterns.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g' )
    END_VERSIONS
    """
}
