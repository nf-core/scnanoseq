//
// Rum UMI Deduplication and optionally split the bam for better parallel processing
//

//
// MODULES
//
include { BAMTOOLS_SPLIT                          } from '../../modules/nf-core/bamtools/split/main'
include { UMITOOLS_DEDUP                          } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP  } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                          } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SPLIT    } from '../../modules/nf-core/samtools/view/main'

include { GROUP_TRANSCRIPTS } from '../../modules/local/group_transcripts'

//
// SUBWORKFLOWS
//
include { BAM_STATS_SAMTOOLS      } from '../../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow DEDUP_UMIS {
    take:
        fasta           // channel: [ val(meta), path(fasta) ]
        fai             // channel: [ val(meta), path(fai) ]
        gtf             // channel: [ val(meta), path(gtf) ]
        in_bam          // channel: [ val(meta), path(bam) ]
        in_bai          // channel: [ val(meta), path(bai) ]
        split_bam       // bool: Split the bam
        genome_aligned  // bool: If the bam is aligned to the genome or not
        dedup_tool      // str: Name of deduplication tool to use
        fasta_delimiter // str: Delimiter character used in the sequence id in fasta

    main:
        ch_versions = Channel.empty()

        ch_undedup_bam = Channel.empty()
        ch_undedup_bai = Channel.empty()

        if (split_bam) {
            ch_split_bam = Channel.empty()

            if (genome_aligned) {
                //
                // MODULE: Bamtools split
                //
                BAMTOOLS_SPLIT ( in_bam )
                ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions.first())
                ch_split_bam = BAMTOOLS_SPLIT.out.bam
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

            } else {
                //
                // MODULE: Group Transcripts
                //
                GROUP_TRANSCRIPTS (
                    fasta,
                    gtf,
                    fasta_delimiter
                )
                ch_versions = ch_versions.mix(GROUP_TRANSCRIPTS.out.versions)

                //
                // MODULE: Samtools View
                //
                SAMTOOLS_VIEW_SPLIT(
                    in_bam
                        .join(in_bai)
                        .combine(GROUP_TRANSCRIPTS.out.grouped_transcripts.flatten())
                        .map{
                            meta, bam, bai, region ->
                                region_basename = region.toString().split('/')[-1]
                                split_region_basename = region_basename.split(/\./)
                                [['id': meta.id + "." + split_region_basename[0]], bam, bai, region]
                        },
                    [[],[]],
                    []
                )
                ch_split_bam = SAMTOOLS_VIEW_SPLIT.out.bam
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_SPLIT.out.versions)
            }

            //
            // SUBWORKFLOW : Sort and index bam
            //
            BAM_SORT_STATS_SAMTOOLS(
                ch_split_bam,
                fasta
            )
            ch_undedup_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_undedup_bai = BAM_SORT_STATS_SAMTOOLS.out.bai
            ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

        }
        else {
            ch_undedup_bam = in_bam
            ch_undedup_bai = in_bai
        }

        ch_dedup_bam = Channel.empty()
        ch_dedup_bai = Channel.empty()

        if (dedup_tool == 'umitools'){
            //
            // MODULE: Umitools Dedup
            //
            UMITOOLS_DEDUP (
                ch_undedup_bam.join(ch_undedup_bai, by: [0]),
                true )
            ch_dedup_bam = UMITOOLS_DEDUP.out.bam
            ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        } else {
            //
            // MODULE: Picard Mark Duplicates
            //
            PICARD_MARKDUPLICATES (
                ch_undedup_bam,
                fasta,
                fai
            )
            ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
            ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
        }

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
                    ch_dedup_bam
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
            ch_dedup_bam = SAMTOOLS_MERGE.out.bam
            ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

            //
            // MODULE: Samtools Index
            //
            SAMTOOLS_INDEX_MERGED( ch_dedup_bam )
            ch_dedup_bai = SAMTOOLS_INDEX_MERGED.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MERGED.out.versions)
        }

        //
        // SUBWORKFLOW: BAM_STATS_SAMTOOLS
        //
        BAM_STATS_SAMTOOLS (
            ch_dedup_bam.join(ch_dedup_bai),
            fasta
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        versions       = ch_versions
        dedup_bam      = ch_dedup_bam
        dedup_bai      = ch_dedup_bai
        dedup_flagstat = BAM_STATS_SAMTOOLS.out.flagstat

        // TODO: Do we need these?
        dedup_stats    = BAM_STATS_SAMTOOLS.out.stats
        dedup_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
}
