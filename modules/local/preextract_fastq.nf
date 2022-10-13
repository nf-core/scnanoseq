process PREEXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    // TODO: We will need a mulled container for this
    conda ( "conda-forge::regex=2022.1.18 conda-forge::biopython=1.79" )

    input:
    tuple val(meta), path(reads)
    path regex_pattern

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
    BC_PATTERN=\$(grep REGEX ${regex_pattern} | sed 's/REGEX: //g')
    FILE_PREFIX=${prefix}

    if [ ${params.split_amount} -gt 0 ]; then
        IDX=\$(basename ${reads} | cut -f2 -d'.')
        FILE_PREFIX=\${FILE_PREFIX}.\${IDX}
    fi

    pre_extract_barcodes.py -i ${reads} \\
                            -r \${BC_PATTERN} \\
                            -o \${FILE_PREFIX}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preextractfastq: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
