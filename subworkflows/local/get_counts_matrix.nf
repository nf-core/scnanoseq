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
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../nf-core/bam_stats_samtools/main'

workflow GET_COUNTS_MATRIX {
    take:
    ch_bam
    ch_gtf

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Count Features
    //
    subread_version = Channel.empty()

    SUBREAD_FEATURECOUNTS ( ch_bam, ch_gtf )
    ch_counts = SUBREAD_FEATURECOUNTS.out.counts
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    //
    // MODULE: Tag Features
    //

    ch_bam_with_counts = Channel.empty()
    ch_bam
        .map{meta, bam -> [meta.id, meta, bam]}
        .join(ch_counts
            .map{meta, bam -> [meta.id, meta, bam]})
        .map{sample_id, bam_meta, bam, counts_meta, counts -> [bam_meta, bam, counts]}
        .set { ch_bam_with_counts }

    TAG_FEATURES ( ch_bam_with_counts )
    ch_tag_bam = TAG_FEATURES.out.feature_bam
    ch_versions = ch_versions.mix(TAG_FEATURES.out.versions)

    //
    // MODULE: Index feature tagged bam
    //
    SAMTOOLS_INDEX ( ch_tag_bam )
    ch_tag_bai = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_tag_bam_bai = ch_tag_bam.join(ch_tag_bai, by: 0)

    //
    // SUBWORKFLOW: Samtools stats
    //

    BAM_STATS_SAMTOOLS ( ch_tag_bam_bai, [] )
    ch_tag_bam_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    //
    // MODULE: Generate the counts matrix
    //

    UMITOOLS_COUNT ( ch_tag_bam_bai )
    ch_count_mtx = UMITOOLS_COUNT.out.counts_matrix
    ch_versions = ch_versions.mix(UMITOOLS_COUNT.out.versions)

    //
    // MODULE: Split the count matrix
    //

    SPLIT_FILE_BY_COLUMN ( ch_count_mtx, 10 )
    ch_split_files = SPLIT_FILE_BY_COLUMN.out.split_files.transpose()
    ch_versions = ch_versions.mix(SPLIT_FILE_BY_COLUMN.out.versions)

    //
    // MODULE: Correct the counts matrix
    //

    CORRECT_COUNTS_MATRIX ( ch_split_files )
    ch_corrected_counts_matrix = CORRECT_COUNTS_MATRIX.out.corrected_counts_matrix
    ch_versions = ch_versions.mix(CORRECT_COUNTS_MATRIX.out.versions)

    //
    // MODULE: Concatenate the matrices column-wise
    //
    MERGE_FILE_BY_COLUMN ( ch_corrected_counts_matrix.groupTuple() )
    ch_counts_mtx = MERGE_FILE_BY_COLUMN.out.col_merged_file
    ch_versions = ch_versions.mix(MERGE_FILE_BY_COLUMN.out.versions)

    emit:
    tag_bam_flagstat = ch_tag_bam_flagstat
    counts_mtx = ch_counts_mtx
    versions = ch_versions

}
