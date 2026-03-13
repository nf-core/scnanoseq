process CORRECT_BARCODES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::editdistance=0.6.0 bioconda::pysam=0.19.1 conda-forge::pygtrie=2.5.0 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' :
        'biocontainers/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' }"

    input:
    tuple val(meta), path(bc_info), path(whitelist), path(bc_count_file)

    output:
    tuple val(meta), path("*.corrected_bc_umi.tsv"), emit: corrected_bc_info
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    correct_barcodes.py \\
        ${args} \\
        --infile ${bc_info} \\
        --outfile ${prefix}.corrected_bc_umi.tsv \\
        --whitelist ${whitelist} \\
        --barcode_count ${bc_count_file} \\
        --skip_header

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.corrected_bc_umi.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
