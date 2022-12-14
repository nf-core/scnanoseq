process PREEXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    // TODO: We will need a mulled container for this
    conda ( "conda-forge::regex=2022.1.18 conda-forge::biopython=1.79" )

    input:
    tuple val(meta), path(reads)
    val regex_pattern

    output:
    tuple val(meta), path("*.R1.fastq"), emit: r1_reads
    tuple val(meta), path("*.R2.fastq"), emit: r2_reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    pre_extract_barcodes.py -i ${reads} \\
                            -r "${regex_pattern}" \\
                            -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preextractfastq: 1.0.0 
    END_VERSIONS
    """
}
