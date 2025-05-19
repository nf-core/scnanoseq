process SPLIT_GTF {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(gtf)

    output:
    path "*.split.gtf"  , emit: split_gtf
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    grep -v '^#' $gtf | awk -F \$'\\t' '{print > \$1".split.gtf"}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch test.split.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
