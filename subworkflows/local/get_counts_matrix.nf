//
//  
//

workflow GET_COUNTS_MATRIX {
    take:
        ch_gtf
        ch_bam

    main:

        // count_featurecounts

        // tag_counts

        ch_counts_mtx = Channel.empty()

    emit:
        ch_counts_mtx
}
