process KALLISTO_BUS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kallisto=0.52.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.52.0--h13ff97a_0' :
        'biocontainers/kallisto:0.52.0--h13ff97a_0' }"

    input:
    tuple val(meta), path(bc_umi_fastq), path(cdna_fastq)
    tuple val(meta2), path(index)
    val technology

    output:
    tuple val(meta), path("bus/output.bus")     , emit: bus
    tuple val(meta), path("bus/matrix.ec")      , emit: ecmap
    tuple val(meta), path("bus/transcripts.txt"), emit: txnames
    tuple val(meta), path("bus/run_info.json")  , emit: run_info
    path "versions.yml"                         , emit: versions_kallisto_bus, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto bus \\
        --long \\
        -i ${index} \\
        -o bus \\
        -x '${technology}' \\
        -t ${task.cpus} \\
        ${args} \\
        ${bc_umi_fastq} \\
        ${cdna_fastq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bus
    touch bus/output.bus bus/matrix.ec bus/transcripts.txt bus/run_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """
}
