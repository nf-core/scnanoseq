//
// Creates the gene-level matrices
//

// Local modules
include { SUBREAD_FEATURECOUNTS } from '../../modules/local/subread_featurecounts'
include { TAG_FEATURES          } from '../../modules/local/tag_features'
include { UMITOOLS_COUNT        } from '../../modules/local/umi_tools_count'
include { SPLIT_FILE_BY_COLUMN  } from '../../modules/local/split_file_by_column'
include { CORRECT_COUNTS_MATRIX } from '../../modules/local/correct_counts_matrix'
include { MERGE_FILE_BY_COLUMN  } from '../../modules/local/merge_file_by_column'

// nf-core modules
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow GET_COUNTS_MATRIX {
    take:
    ch_bam
    ch_gtf

    main:

    //
    // MODULE: Count Features 
    //
    subread_version = Channel.empty()

    SUBREAD_FEATURECOUNTS ( ch_bam, ch_gtf )
    ch_counts = SUBREAD_FEATURECOUNTS.out.counts

    subread_version = SUBREAD_FEATURECOUNTS.out.versions

    //
    // MODULE: Tag Features
    //
    TAG_FEATURES ( ch_bam.join(ch_counts, by: 0) )
    ch_tag_bam = TAG_FEATURES.out.feature_bam

    //
    // MODULE: Index feature tagged bam
    //
    SAMTOOLS_INDEX (ch_tag_bam)
    ch_tag_bam_bai = SAMTOOLS_INDEX.out.bai

    //
    // MODULE: Generate the counts matrix
    //

    UMITOOLS_COUNT ( ch_tag_bam.join(ch_tag_bam_bai, by: 0) )
    ch_count_mtx = UMITOOLS_COUNT.out.counts_matrix

    //
    // MODULE: Split the count matrix
    //

    SPLIT_FILE_BY_COLUMN ( ch_count_mtx, 10 )
    ch_split_files = SPLIT_FILE_BY_COLUMN.out.split_files.transpose()

    //
    // MODULE: Correct the counts matrix
    //

    CORRECT_COUNTS_MATRIX ( ch_split_files )
    ch_corrected_counts_matrix = CORRECT_COUNTS_MATRIX.out.corrected_counts_matrix

    //
    // MODULE: Concatenate the matrices column-wise
    //
    MERGE_FILE_BY_COLUMN ( ch_corrected_counts_matrix.groupTuple() )
    counts_mtx = MERGE_FILE_BY_COLUMN.out.col_merged_file

    emit:
    counts_mtx
    subread_version

}
