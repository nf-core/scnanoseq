process KALLISTO_INDEX {
    tag "kallisto_index"
    label 'process_high'

    conda "bioconda::kallisto=0.52.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.52.0--h13ff97a_0' :
        'biocontainers/kallisto:0.52.0--h13ff97a_0' }"

    input:
    tuple val(meta), path(transcript_fasta), path(genome_fasta)

    output:
    tuple val(meta), path("*.idx"), emit: index
    path "versions.yml"           , emit: versions_kallisto_index, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "kallisto"
    // The genome is used as a d-list so that reads originating outside the
    // transcriptome are discarded rather than misassigned.
    def d_list = genome_fasta ? "-d ${genome_fasta}" : ''
    """
    kallisto index \\
        -i ${prefix}.idx \\
        -t ${task.cpus} \\
        ${d_list} \\
        ${args} \\
        ${transcript_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "kallisto"
    """
    touch ${prefix}.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed 's/^kallisto, version //')
    END_VERSIONS
    """
}
