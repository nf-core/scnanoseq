process SORT_GTF {
    tag "$meta.id"
    label 'process_low'

//    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
//        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.sorted.gtf"), emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${gtf} | awk '\$1 !~ /^#/ { print \$0 }' | sort -k1,1 -k4,4n -k5,5n -V > ${prefix}.sorted.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
        sort: \$(echo \$(sort --version) | sed 's/^.*sort (GNU coreutils) //; s/ .*//')
    END_VERSIONS

    """
}
