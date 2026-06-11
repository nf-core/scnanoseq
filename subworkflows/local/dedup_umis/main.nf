//
// Rum UMI Deduplication and optionally split the bam for better parallel processing
//

//
// MODULES
//
include { BAMTOOLS_SPLIT                          } from '../../../modules/nf-core/bamtools/split/main'
include { UMITOOLS_DEDUP                          } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                          } from '../../../modules/nf-core/samtools/merge/main'
include { SPLIT_BAM                               } from '../../../modules/local/split_bam'
include { GROUP_TRANSCRIPTS                       } from '../../../modules/local/group_transcripts'
include { PICARD_MARKDUPLICATES                   } from '../../../modules/nf-core/picard/markduplicates/main'

//
// SUBWORKFLOWS
//
include { BAM_STATS_SAMTOOLS      } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow DEDUP_UMIS {
    take:
        ch_fasta // channel: [ val(meta), path(fasta) ]
        ch_fai   // channel: [ val(meta), path(fai) ]
        ch_gtf   // channel: [ val(meta), path(gtf) ]
        ch_bam   // channel: [ val(meta), path(bam) ]
        ch_bai   // channel: [ val(meta), path(bai) ]

        val_split_bam       // bool: Split the bam
        val_genome_aligned  // bool: If the bam is aligned to the genome or not
        val_dedup_tool      // str: Name of deduplication tool to use
        val_fasta_delimiter // str: Delimiter character used in the sequence id in fasta

    main:
        ch_versions = channel.empty()

        ch_undedup_bam = channel.empty()
        ch_undedup_bai = channel.empty()

        if (val_split_bam) {
            ch_split_bam = channel.empty()

            if (val_genome_aligned) {
                //
                // MODULE: Bamtools split
                //
                BAMTOOLS_SPLIT ( ch_bam )
                ch_split_bam = BAMTOOLS_SPLIT.out.bam
                    .map{
                        _meta, bam ->
                            [bam]
                    }
                    .flatten()
                    .map{
                        bam ->
                            def bam_basename = bam.toString().split('/')[-1]
                            def split_bam_basename = bam_basename.split(/\./)
                            def new_meta = [
                                'id': split_bam_basename.take(split_bam_basename.size()-1).join("."),
                            ]
                            [ new_meta, bam ]
                    }

            } else {
                //
                // MODULE: Group Transcripts
                //
                GROUP_TRANSCRIPTS (
                    ch_fasta,
                    ch_gtf,
                    val_fasta_delimiter
                )
                ch_versions = ch_versions.mix(GROUP_TRANSCRIPTS.out.versions_group_transcripts)

                //
                // MODULE: Samtools View
                //
                SPLIT_BAM(
                    ch_bam
                        .join(ch_bai)
                        .combine(GROUP_TRANSCRIPTS.out.grouped_transcripts.flatten())
                        .map{
                            meta, bam, bai, region ->
                                def region_basename = region.toString().split('/')[-1]
                                def split_region_basename = region_basename.split(/\./)
                                [['id': meta.id + "." + split_region_basename[0]], bam, bai, region]
                        }
                )
                ch_split_bam = SPLIT_BAM.out.split_bam
                ch_versions = ch_versions.mix(SPLIT_BAM.out.versions_split_bam)
            }

            //
            // SUBWORKFLOW : Sort and index bam
            //
            BAM_SORT_STATS_SAMTOOLS(
                ch_split_bam,
                ch_fasta.first()
            )

            ch_undedup_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_undedup_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

        }
        else {
            ch_undedup_bam = ch_bam
            ch_undedup_bai = ch_bai
        }

        ch_dedup_bam = channel.empty()
        ch_dedup_bai = channel.empty()

        if (val_dedup_tool == 'umitools'){
            //
            // MODULE: Umitools Dedup
            //
            UMITOOLS_DEDUP (
                ch_undedup_bam.join(ch_undedup_bai, by: [0]),
                true )
            ch_dedup_bam = UMITOOLS_DEDUP.out.bam

        } else {
            //
            // MODULE: Picard Mark Duplicates
            //
            PICARD_MARKDUPLICATES (
                ch_undedup_bam,
                ch_fasta.first(),
                ch_fai.first()
            )
            ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        }

        //
        // MODULE: Samtools Index
        //
        SAMTOOLS_INDEX_DEDUP( UMITOOLS_DEDUP.out.bam )
        ch_dedup_bai = SAMTOOLS_INDEX_DEDUP.out.bai

        if (val_split_bam) {
            //
            // MODULE: Samtools Merge
            //
            SAMTOOLS_MERGE (
                    ch_dedup_bam
                        .map{
                            _meta, bam ->
                                def bam_basename = bam.toString().split('/')[-1]
                                def split_bam_basename = bam_basename.split(/\./)
                                def new_meta = [ 'id': split_bam_basename[0] ]
                            [ new_meta, bam ]
                        }
                        .groupTuple()
                        .join(
                            ch_dedup_bai
                                .map{
                                    _meta, bai ->
                                        def bai_basename = bai.toString().split('/')[-1]
                                        def split_bai_basename = bai_basename.split(/\./)
                                        def new_meta = [ 'id': split_bai_basename[0] ]
                                    [ new_meta, bai ]
                                }
                                .groupTuple()
                        ),
                ch_fasta
                    .join(ch_fai)
                    .map { meta, fasta_file, fai_file ->
                        [meta, fasta_file, fai_file, "$projectDir/assets/dummy_file.txt"]
                    }
                    .first()

            )
            ch_dedup_bam = SAMTOOLS_MERGE.out.bam

            //
            // MODULE: Samtools Index
            //
            SAMTOOLS_INDEX_MERGED( ch_dedup_bam )
            ch_dedup_bai = SAMTOOLS_INDEX_MERGED.out.bai
        }

        //
        // SUBWORKFLOW: BAM_STATS_SAMTOOLS
        //
        BAM_STATS_SAMTOOLS (
            ch_dedup_bam.join(ch_dedup_bai),
            ch_fasta.first()
        )

    emit:
        versions       = ch_versions                     // channel: [ val(meta), path(versions) ]
        dedup_bam      = ch_dedup_bam                    // channel: [ val(meta), path(bam) ]
        dedup_bai      = ch_dedup_bai                    // channel: [ val(meta), path(bai) ]
        dedup_flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
}
