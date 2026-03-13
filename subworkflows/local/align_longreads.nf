//
// Performs alignment
//

// SUBWORKFLOWS
include { BAM_SORT_STATS_SAMTOOLS                                     } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FILTERED } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'

// MODULES
include { MINIMAP2_INDEX                          } from '../../modules/nf-core/minimap2/index'
include { MINIMAP2_ALIGN                          } from '../../modules/nf-core/minimap2/align'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_MAPPED } from '../../modules/nf-core/samtools/view'

include { RSEQC_READDISTRIBUTION } from '../../modules/nf-core/rseqc/readdistribution/main'
include { NANOCOMP               } from '../../modules/nf-core/nanocomp/main'


workflow ALIGN_LONGREADS {
    take:
        fasta       // channel: [ val(meta), path(fasta) ]
        fai         // channel: [ val(meta), path(fai) ]
        gtf         // channel: [ val(meta), path(gtf) ]
        fastq       // channel: [ val(meta), path(fastq) ]
        rseqc_bed   // channel: [ val(meta), path(rseqc_bed) ]

        skip_save_minimap2_index // bool: Skip saving the minimap2 index
        skip_qc                  // bool: Skip qc steps
        skip_rseqc               // bool: Skip RSeQC
        skip_bam_nanocomp        // bool: Skip Nanocomp

    main:
        ch_versions = Channel.empty()
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
        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        BAM_SORT_STATS_SAMTOOLS ( MINIMAP2_ALIGN.out.bam, fasta )
        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

        // acquire only mapped reads from bam for downstream processing
        // NOTE: some QCs steps are performed on the full BAM
        SAMTOOLS_FILTER_MAPPED (
            BAM_SORT_STATS_SAMTOOLS.out.bam
                .join( BAM_SORT_STATS_SAMTOOLS.out.bai, by: 0 )
                .combine(["$projectDir/assets/dummy_file.txt"]),
            [[],[]],
            []
        )

        ch_minimap_mapped_only_bam = SAMTOOLS_FILTER_MAPPED.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_FILTER_MAPPED.out.versions)

        BAM_SORT_STATS_SAMTOOLS_FILTERED (
            ch_minimap_mapped_only_bam,
            fasta
        )

        ch_minimap_filtered_sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam
        ch_minimap_filtered_sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai
        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_FILTERED.out.versions)

        //
        // MODULE: RSeQC read distribution for BAM files (unfiltered for QC purposes)
        //
        ch_rseqc_read_dist = Channel.empty()
        if (!skip_qc && !skip_rseqc) {
            RSEQC_READDISTRIBUTION ( BAM_SORT_STATS_SAMTOOLS.out.bam, rseqc_bed )
            ch_rseqc_read_dist = RSEQC_READDISTRIBUTION.out.txt
            ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions)
        }

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
        versions = ch_versions

        // Bam and Bai
        sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam
        sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai

        // SAMtool stats from initial mapping
        stats = BAM_SORT_STATS_SAMTOOLS.out.stats
        flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats

        // RSeQC stats
        rseqc_read_dist = ch_rseqc_read_dist

        // Nanoplot stats
        nanocomp_bam_html = ch_nanocomp_bam_html
        nanocomp_bam_txt = ch_nanocomp_bam_txt
}
