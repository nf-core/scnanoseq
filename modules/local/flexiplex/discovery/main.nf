process FLEXIPLEX_DISCOVERY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiplex:1.01--py310h84f13bb_1':
        'biocontainers/flexiplex:1.01--py310h84f13bb_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("flexiplex_barcodes_counts.txt")  , emit: barcode_counts
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}${meta.part ? "_part_${meta.part}" : ''}"
    """
    # Run in discovery mode
    flexiplex \\
        ${args} \\
        -p ${task.cpus} \\
        ${reads}  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch flexiplex_barcodes_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """
}
