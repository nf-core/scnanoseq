process PREEXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_low'
    //stageInMode 'copy'

    conda "conda-forge::regex=2022.1.18 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-497e12343495a3e0f3b3459542cc8ad37575d9fa:483e027ac6835fcb80b9cfef4de8c89b67343941-0' :
        'biocontainers/mulled-v2-497e12343495a3e0f3b3459542cc8ad37575d9fa:483e027ac6835fcb80b9cfef4de8c89b67343941-0' }"

    input:
    tuple val(meta), path(reads), path(bc_list)
    val bc_format

    output:
    tuple val(meta), path("*.putative_bc_umi.tsv"), emit: barcode_info
    tuple val(meta), path("*.extracted.fastq"), emit: extracted_fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pre_extract_barcodes.py \\
        -i ${reads} \\
        -b ${bc_list} \\
        -o ${prefix}.extracted \\
        -f ${bc_format} \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.putative_bc_umi.tsv
    touch ${prefix}.extracted.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
