process FLEXIPLEX_ASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiplex:1.01--py310h84f13bb_1':
        'biocontainers/flexiplex:1.01--py310h84f13bb_1' }"

    input:
    tuple val(meta), path(reads), path(barcodes)

    output:
    tuple val(meta), path("*flexiplex.fastq")               , emit: reads
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}${meta.part ? "_part_${meta.part}" : ''}"
    """
    # Run in assignment mode

    flexiplex \\
        ${args} \\
        -k ${barcodes} \\
        -p ${task.cpus} \\
        ${reads} \\
        > ${prefix}.flexiplex.fastq
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.flexiplex.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """
}
