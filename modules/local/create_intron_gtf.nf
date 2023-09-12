process CREATE_INTRON_GTF {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.intron.gtf"), emit: intron_gtf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat ${gtf} | \\
        tr -d '\\000'| \\
        awk -F \$'\\t' 'BEGIN{OFS="\\t"} \$3="intron", \$7="+"'  | \\
        grep -v exon_id > ${prefix}.intron.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version) | sed 's/^.*GNU Awk //; s/ .*//')
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
        grep: \$(echo \$(grep --version) | sed 's/^.*grep (GNU grep) //; s/ .*//')
        tr: \$(echo \$(tr --version) | sed 's/^.*tr (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
