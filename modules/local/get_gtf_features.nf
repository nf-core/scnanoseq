process GET_GTF_FEATURES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(gtf)
    val feature

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk 'OFS="\\t", \$1 !~ /^#/ {if (\$3 == "exon") print \$0}' ${gtf} > ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version) | sed 's/^.*GNU Awk //; s/ .*//')
    END_VERSIONS
    """
}
