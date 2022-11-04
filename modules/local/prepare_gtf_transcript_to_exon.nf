process TRANSCRIPT_TO_EXON {
    tag "$gtf"
    label 'process_low'

    //TODO: Figure out what environment is needed

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
