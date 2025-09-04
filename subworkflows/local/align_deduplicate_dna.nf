//
// Performs alignment and deduplication for DNA samples
//

// MODULES
include { MINIMAP2_INDEX                          } from '../../modules/nf-core/minimap2/index'
include { MINIMAP2_ALIGN                          } from '../../modules/nf-core/minimap2/align'
include { PICARD_MARKDUPLICATES                   } from '../../modules/nf-core/picard/markduplicates'
include { BAM_SORT_STATS_SAMTOOLS                 } from '../../subworkflows/nf-core/bam_sort_stats_samtools'
include { FLEXIFORMATTER                          } from '../../modules/local/flexiformatter'
include { NANOCOMP                                } from '../../modules/nf-core/nanocomp/main'

workflow ALIGN_DEDUPLICATE_DNA {
    take:
        fasta           // channel: [ val(meta), path(fasta) ]
        fai             // channel: [ val(meta), path(fai) ]
        fastq           // channel: [ val(meta), path(fastq) ]

        skip_save_minimap2_index // bool: Skip saving the minimap2 index
        skip_qc                  // bool: Skip qc steps
        skip_bam_nanocomp        // bool: Skip Nanocomp
        skip_dedup               // bool: Skip deduplication

    main:
        ch_versions              = Channel.empty()

        // Minimap results
        minimap_bam              = Channel.empty()
        minimap_bai              = Channel.empty()

        // Deduplicated bam file
        dedup_bam                = Channel.empty()
        dedup_bai                = Channel.empty()

        // SAMtool stats after dedup
        stats                    = Channel.empty()
        flagstat                 = Channel.empty()
        idxstats                 = Channel.empty()

        // NanoComp results
        nanocomp_bam_html        = Channel.empty()
        nanocomp_bam_txt         = Channel.empty()

        //
        // MINIMAP2_INDEX
        //
        if (skip_save_minimap2_index) {
            MINIMAP2_INDEX ( fasta )
            ch_minimap_ref = MINIMAP2_INDEX.out.index
            ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        } else {
            ch_minimap_ref = fasta
        }

        //
        // MINIMAP2_ALIGN
        //

        MINIMAP2_ALIGN (
            fastq,
            ch_minimap_ref,
            true,
            "bai",
            "",
            ""
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        //
        // MODULE: MarkDuplicates
        //
        final_bam = MINIMAP2_ALIGN.out.bam
        if( !skip_dedup ) {
            //
            // MODULE: Picard Mark Duplicates
            //
            PICARD_MARKDUPLICATES (
                MINIMAP2_ALIGN.out.bam,
                fasta,
                fai
            )
            final_bam = PICARD_MARKDUPLICATES.out.bam
            ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
        }


        //
        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        // TODO: No reason that this is again sorting and indexing.
        // Change to STATS_SAMTOOLS
        BAM_SORT_STATS_SAMTOOLS ( final_bam, fasta )
        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

        //
        // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
        //
        ch_nanocomp_bam_html = Channel.empty()
        ch_nanocomp_bam_txt = Channel.empty()

        if (!skip_qc && !skip_bam_nanocomp) {

            NANOCOMP (
                BAM_SORT_STATS_SAMTOOLS.out.bam
                    .collect{it[1]}
                    .map{
                        [ [ 'id': 'nanocomp_bam.' ] , it ]
                    }
            )

            ch_nanocomp_bam_html = NANOCOMP.out.report_html
            ch_nanocomp_bam_txt = NANOCOMP.out.stats_txt
            ch_versions = ch_versions.mix( NANOCOMP.out.versions )
        }

    emit:
        // Versions
        versions                 = ch_versions

        // Minimap results
        minimap_bam              = MINIMAP2_ALIGN.out.bam
        minimap_bai              = MINIMAP2_ALIGN.out.index

        // Deduplicated bam file
        dedup_bam                = BAM_SORT_STATS_SAMTOOLS.out.bam
        dedup_bai                = BAM_SORT_STATS_SAMTOOLS.out.bai

        // SAMtool stats after dedup
        stats                    = BAM_SORT_STATS_SAMTOOLS.out.stats
        flagstat                 = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        idxstats                 = BAM_SORT_STATS_SAMTOOLS.out.idxstats

        // NanoComp results
        nanocomp_bam_html        = ch_nanocomp_bam_html
        nanocomp_bam_txt         = ch_nanocomp_bam_txt
}
