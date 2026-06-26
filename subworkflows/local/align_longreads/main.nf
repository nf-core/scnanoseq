//
// Performs alignment
//

// SUBWORKFLOWS
include { BAM_SORT_STATS_SAMTOOLS                                     } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FILTERED } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

// MODULES
include { MINIMAP2_INDEX                          } from '../../../modules/nf-core/minimap2/index'
include { MINIMAP2_ALIGN                          } from '../../../modules/nf-core/minimap2/align'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_MAPPED } from '../../../modules/nf-core/samtools/view'

include { RSEQC_READDISTRIBUTION } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { NANOCOMP               } from '../../../modules/nf-core/nanocomp/main'


workflow ALIGN_LONGREADS {
    take:
        ch_fasta     // channel: [ val(meta), path(fasta) ]
        ch_fastq     // channel: [ val(meta), path(fastq) ]
        ch_rseqc_bed // channel: [ val(meta), path(rseqc_bed) ]

        val_skip_save_minimap2_index // bool: Skip saving the minimap2 index
        val_skip_qc                  // bool: Skip qc steps
        val_skip_rseqc               // bool: Skip RSeQC
        val_skip_bam_nanocomp        // bool: Skip Nanocomp

    main:
        //
        // MINIMAP2_INDEX
        //
        if (val_skip_save_minimap2_index) {
            MINIMAP2_INDEX ( ch_fasta )
            ch_minimap_ref = MINIMAP2_INDEX.out.index
        } else {
            ch_minimap_ref = ch_fasta
        }

        //
        // MINIMAP2_ALIGN
        //
        MINIMAP2_ALIGN (
            ch_fastq,
            ch_minimap_ref.first(),
            true,
            "bai",
            "",
            ""
        )

        //
        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        BAM_SORT_STATS_SAMTOOLS ( MINIMAP2_ALIGN.out.bam, ch_fasta.first() )

        // acquire only mapped reads from bam for downstream processing
        // NOTE: some QCs steps are performed on the full BAM
        SAMTOOLS_FILTER_MAPPED (
            BAM_SORT_STATS_SAMTOOLS.out.bam
                .join( BAM_SORT_STATS_SAMTOOLS.out.bai, by: 0 )
                .combine(["$projectDir/assets/dummy_file.txt"]),
            [[],[], []],
            [],
            []
        )

        ch_minimap_mapped_only_bam = SAMTOOLS_FILTER_MAPPED.out.bam

        BAM_SORT_STATS_SAMTOOLS_FILTERED (
            ch_minimap_mapped_only_bam,
            ch_fasta.first()
        )

        //
        // MODULE: RSeQC read distribution for BAM files (unfiltered for QC purposes)
        //
        ch_rseqc_read_dist = channel.empty()
        if (!val_skip_qc && !val_skip_rseqc) {
            RSEQC_READDISTRIBUTION (
                BAM_SORT_STATS_SAMTOOLS.out.bam.join( BAM_SORT_STATS_SAMTOOLS.out.bai, by: 0 ),
                ch_rseqc_bed
            )
            ch_rseqc_read_dist = RSEQC_READDISTRIBUTION.out.txt
        }

        //
        // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
        //
        ch_nanocomp_bam_html = channel.empty()
        ch_nanocomp_bam_txt = channel.empty()

        if (!val_skip_qc && !val_skip_bam_nanocomp) {

            NANOCOMP (
                BAM_SORT_STATS_SAMTOOLS.out.bam
                    .collect{ v -> v[1] }
                    .map{ bams ->
                        [ [ 'id': 'nanocomp_bam.' ] , bams ]
                    }
            )

            ch_nanocomp_bam_html = NANOCOMP.out.report_html
            ch_nanocomp_bam_txt = NANOCOMP.out.stats_txt
        }

    emit:
        // Bam and Bai
        sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam // channel: [ val(meta), path(bam) ]
        sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai // channel: [ val(meta), path(bai) ]

        // SAMtool stats from initial mapping
        stats = BAM_SORT_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), path(stats) ]
        flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
        idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

        // RSeQC stats
        rseqc_read_dist = ch_rseqc_read_dist // channel: [ val(meta), path(rseqc_read_dist) ]

        // Nanoplot stats
        nanocomp_bam_html = ch_nanocomp_bam_html // channel: [ val(meta), path(nanocomp_bam_html) ]
        nanocomp_bam_txt = ch_nanocomp_bam_txt   // channel: [ val(meta), path(nanocomp_bam_txt) ]
}
