//
// Rum UMI Dedupliation and optionally split the bam for better parallel processing
//

//
// MODULES
//
include { BAMTOOLS_SPLIT                          } from '../../modules/nf-core/bamtools/split/main'
include { UMITOOLS_DEDUP                          } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SPLIT  } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP  } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                          } from '../../modules/nf-core/samtools/merge/main'

//
// SUBWORKFLOWS
//
include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow UMITOOLS_DEDUP_SPLIT {
    take:
        fasta       // channel: [ val(meta), path(fasta) ]
        fai         // channel: [ val(meta), path(fai) ]
        in_bam      // channel: [ val(meta), path(bam) ]
        in_bai      // channel: [ val(meta), path(bai) ]
        split_bam   // bool: Split the bam

    main:
        ch_versions = Channel.empty()

        if (split_bam) {
            //
            // MODULE: Bamtools Split
            //
            BAMTOOLS_SPLIT ( in_bam )
            ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions.first())
            ch_undedup_bam = BAMTOOLS_SPLIT.out.bam
                .map{
                    meta, bam ->
                        [bam]
                }
                .flatten()
                .map{
                    bam ->
                        bam_basename = bam.toString().split('/')[-1]
                        split_bam_basename = bam_basename.split(/\./)
                        meta = [
                            'id': split_bam_basename.take(split_bam_basename.size()-1).join("."),
                        ]
                        [ meta, bam ]
                }
            //
            // MODULE: Samtools Index
            //
            SAMTOOLS_INDEX_SPLIT( ch_undedup_bam )
            ch_undedup_bai = SAMTOOLS_INDEX_SPLIT.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SPLIT.out.versions.first())
        }
        else {
            ch_undedup_bam = in_bam
            ch_undedup_bai = in_bai
        }

        //
        // MODULE: Umitools Dedup
        //
        UMITOOLS_DEDUP ( ch_undedup_bam.join(ch_undedup_bai, by: [0]), true )
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        //
        // MODULE: Samtools Index
        //
        SAMTOOLS_INDEX_DEDUP( UMITOOLS_DEDUP.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions)

        if (split_bam) {
            //
            // MODULE: Samtools Merge
            //
            SAMTOOLS_MERGE (
                UMITOOLS_DEDUP.out.bam
                    .map{
                        meta, bam ->
                            bam_basename = bam.toString().split('/')[-1]
                            split_bam_basename = bam_basename.split(/\./)
                            meta = [ 'id': split_bam_basename[0] ]
                        [ meta, bam ]
                    }
                    .groupTuple(),
                fasta,
                fai)
            ch_dedup_single_bam = SAMTOOLS_MERGE.out.bam
            ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

            //
            // MODULE: Samtools Index
            //
            SAMTOOLS_INDEX_MERGED( ch_dedup_single_bam )
            ch_dedup_single_bai = SAMTOOLS_INDEX_MERGED.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MERGED.out.versions)
        }
        else {
            ch_dedup_single_bam = UMITOOLS_DEDUP.out.bam
            ch_dedup_single_bai = SAMTOOLS_INDEX_DEDUP.out.bai
        }

        //
        // SUBWORKFLOW: BAM_STATS_SAMTOOLS
        //
        BAM_STATS_SAMTOOLS (
            ch_dedup_single_bam.join(ch_dedup_single_bai),
            fasta
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        versions       = ch_versions
        dedup_bam      = UMITOOLS_DEDUP.out.bam
        dedup_log      = UMITOOLS_DEDUP.out.log
        dedup_bai      = SAMTOOLS_INDEX_DEDUP.out.bai
        dedup_flagstat = BAM_STATS_SAMTOOLS.out.flagstat

        // TODO: Do we need these?
        dedup_stats    = BAM_STATS_SAMTOOLS.out.stats
        dedup_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
}
