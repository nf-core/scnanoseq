process GROUP_TRANSCRIPTS {
    tag "group_transcripts"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(gtf)
    val delimiter

    output:
    path "*.transcripts.txt", emit: grouped_transcripts
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def delim_arg = delimiter ? "-d '${delimiter}'" : ''

    """
    group_transcripts.sh \\
        -f ${fasta} \\
        -g ${gtf} \\
        ${delim_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch chr1.transcripts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
