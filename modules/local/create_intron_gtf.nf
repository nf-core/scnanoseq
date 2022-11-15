process CREATE_INTRON_GTF {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-genepredtogtf:377--h0b8a92a_4' :
        'quay.io/biocontainers/ucsc-genepredtogtf:377--h0b8a92a_4' }"

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
        create_intron_gtf: 1.0 
    END_VERSIONS
    """
}
