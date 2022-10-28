process TAG_BARCODES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(r1_fastq), val(bc_length), val(umi_length)

    output:
    tuple val(meta), path("*.tagged.bam"), emit: tagged_bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    tag_barcodes.py \\
        -b ${bam} \\
        -f ${r1_fastq} \\
        -o ${prefix}.tagged.bam \\
        -bl ${bc_length} \\
        -ul ${umi_length}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tag_bam: v1.0 
    END_VERSIONS
    """
}
