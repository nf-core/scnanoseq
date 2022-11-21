process TRANSCRIPT_TO_EXON {
    tag "$gtf"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path gtf

    output:
    path("processed.gtf"), emit: ch_processed_gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    awk 'BEGIN{FS="\t"; OFS="\t"} \$3 == "transcript" { \$3="exon"; print}' $gtf > processed.gtf
    """
}
