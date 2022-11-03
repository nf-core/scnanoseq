//
//
//


// Local modules
include { SORT_GTF } from '../../modules/local/sort_gtf'

// nf-core modules
include { STRINGTIE_STRINGTIE } from '../../modules/nf-core/stringtie/stringtie/main'

workflow GET_COUNTS_MATRIX {
    take:
        ch_gtf
        ch_bam

    main:
        // make_gtf_stringtie2
        STRINGTIE_STRINGTIE( ch_bam, ch_gtf )
        ch_stringtie_gtf = STRINGTIE_STRINGTIE.out.gtf

        // sort_gtf
        SORT_GTF( ch_stringtie_gtf )
        ch_sorted_stringtie_gtf = SORT_GTF.out.gtf

        // merge_gtf_stringtie2

        // count_featurecounts

        // tag_counts

        ch_counts_mtx = Channel.empty()

    emit:
        ch_counts_mtx
}
