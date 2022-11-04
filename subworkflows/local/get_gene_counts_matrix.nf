//
// Creates the gene-level matrices
//

// Local modules
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE         } from '../../modules/local/subread_featurecounts'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE_INTRON2 } from '../../modules/local/subread_featurecounts'

workflow GET_GENE_COUNTS_MATRIX {
    take:
    ch_bam
    ch_gtf
    gtf_preparation_method
    ch_gtf_intron2

    main:

    // count_featurecounts
    subread_version = Channel.empty()

    SUBREAD_FEATURECOUNTS_GENE ( ch_bam, ch_gtf )
    gene_counts_mtx = SUBREAD_FEATURECOUNTS_GENE.out.counts

    //TODO: enable this once we have intron 2 GTF in the workflow
    // for intron method 2, perform a second round for introns count
    //if (gtf_preparation_method == "2") {
      //  SUBREAD_FEATURECOUNTS_GENE_INTRON2 ( ch_bam, ch_gtf_intron2 )
        //gene_counts_intron2_mtx = SUBREAD_FEATURECOUNTS_GENE_INTRON2.out.counts
        //}

    subread_version = SUBREAD_FEATURECOUNTS_GENE.out.versions

    // tag_counts

    tag_gene_counts_mtx = Channel.empty()
    tag_gene_counts_intron2_mtx = Channel.empty()

    emit:
    gene_counts_mtx //TODO: can be removed once tagging is added (just here fore test)
    tag_gene_counts_mtx
    tag_gene_counts_intron2_mtx
    subread_version

}
