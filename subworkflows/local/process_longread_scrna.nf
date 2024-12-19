//
// Performs alignment
//

// SUBWORKFLOWS
include { ALIGN_LONGREADS         } from '../../subworkflows/local/align_longreads'
include { QUANTIFY_SCRNA_ISOQUANT } from '../../subworkflows/local/quantify_scrna_isoquant'
include { QUANTIFY_SCRNA_OARFISH  } from '../../subworkflows/local/quantify_scrna_oarfish'
include { UMITOOLS_DEDUP_SPLIT    } from '../../subworkflows/local/umitools_dedup_split'

// MODULES
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TAGGED       } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TAGGED } from '../../modules/nf-core/samtools/flagstat'

include { TAG_BARCODES } from '../../modules/local/tag_barcodes'


workflow PROCESS_LONGREAD_SCRNA {
    take:
        fasta        // channel: [ val(meta), path(fasta) ]
        fai          // channel: [ val(meta), path(fai) ]
        gtf          // channel: [ val(meta), path(gtf) ]
        fastq        // channel: [ val(meta), path(fastq) ]
        rseqc_bed    // channel: [ val(meta), path(rseqc_bed) ]
        read_bc_info // channel: [ val(meta), path(read_barcode_info) ]

        skip_save_minimap2_index // bool: Skip saving the minimap2 index
        skip_qc                  // bool: Skip qc steps
        skip_rseqc               // bool: Skip RSeQC
        skip_bam_nanocomp        // bool: Skip Nanocomp
        quantifier               // str: Quantifier
        skip_seurat              // bool: Skip seurat qc steps
        split_umitools_bam       // bool: Skip splitting on chromsome for umitools

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Align long Read Data
        //

        ALIGN_LONGREADS(
            fasta,
            fai,
            gtf,
            fastq,
            rseqc_bed,
            skip_save_minimap2_index,
            skip_qc,
            skip_rseqc,
            skip_bam_nanocomp
        )
        ch_versions = ch_versions.mix(ALIGN_LONGREADS.out.versions)

        //
        // MODULE: Tag Barcodes
        //

        TAG_BARCODES (
            ALIGN_LONGREADS.out.sorted_bam
                .join( ALIGN_LONGREADS.out.sorted_bai, by: 0 )
                .join( read_bc_info, by: 0)
        )
        ch_versions = ch_versions.mix(TAG_BARCODES.out.versions)

        //
        // MODULE: Index Tagged Bam
        //
        SAMTOOLS_INDEX_TAGGED ( TAG_BARCODES.out.tagged_bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TAGGED.out.versions)

        //
        // MODULE: Flagstat Tagged Bam
        //
        SAMTOOLS_FLAGSTAT_TAGGED (
            TAG_BARCODES.out.tagged_bam
                .join( SAMTOOLS_INDEX_TAGGED.out.bai, by: [0])
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_TAGGED.out.versions)

        //
        // SUBWORKFLOW: UMI Deduplication
        //
        if (!params.skip_dedup && quanitifer.equals("isoquant")) {
            UMITOOLS_DEDUP_SPLIT(
                fasta,
                fai,
                TAG_BARCODES.out.tagged_bam,
                SAMTOOLS_INDEX_TAGGED.out.bai,
                split_umitools_bam
            )

            ch_versions = ch_versions.mix(UMITOOLS_DEDUP_SPLIT.out.versions)
        }

        //
        // SUBWORKFLOW: Quantify Features
        //

        ch_gene_qc_stats = Channel.empty()
        ch_transcript_qc_stats = Channel.empty()

        if (quantifier.equals("oarfish")) {
            QUANTIFY_SCRNA_OARFISH (
                params.skip_dedup? TAG_BARCODES.out.tagged_bam : UMITOOLS_DEDUP_SPLIT.out.dedup_bam,
                params.skip_dedup? SAMTOOLS_INDEX_TAGGED.out.bai : UMITOOLS_DEDUP_SPLIT.out.dedup_bai,
                params.skip_dedup? SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat : UMITOOLS_DEDUP_SPLIT.out.dedup_flagstat,
                fasta,
                skip_qc,
                skip_seurat
            )
            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_OARFISH.out.versions)
            ch_transcript_qc_stats = QUANTIFY_SCRNA_OARFISH.out.transcript_qc_stats
        } else {
            QUANTIFY_SCRNA_ISOQUANT (
                params.skip_dedup? TAG_BARCODES.out.tagged_bam : UMITOOLS_DEDUP_SPLIT.out.dedup_bam,
                params.skip_dedup? SAMTOOLS_INDEX_TAGGED.out.bai : UMITOOLS_DEDUP_SPLIT.out.dedup_bai,
                params.skip_dedup? SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat : UMITOOLS_DEDUP_SPLIT.out.dedup_flagstat,
                fasta,
                fai,
                gtf,
                skip_qc,
                skip_seurat
            )

            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_ISOQUANT.out.versions)
            ch_gene_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.gene_qc_stats
            ch_transcript_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.transcript_qc_stats
        }

    emit:
        // Versions
        versions = ch_versions

        minimap_bam = ALIGN_LONGREADS.out.sorted_bam
        minimap_bai = ALIGN_LONGREADS.out.sorted_bai
        minimap_stats = ALIGN_LONGREADS.out.stats
        minimap_flagstat = ALIGN_LONGREADS.out.flagstat
        minimap_idxstats = ALIGN_LONGREADS.out.idxstats
        minimap_rseqc_read_dist = ALIGN_LONGREADS.out.rseqc_read_dist
        minimap_nanocomp_bam_txt = ALIGN_LONGREADS.out.nanocomp_bam_txt

        bc_tagged_bam = TAG_BARCODES.out.tagged_bam
        bc_tagged_bai = SAMTOOLS_INDEX_TAGGED.out.bai
        bc_tagged_flagstat = SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat

        dedup_bam = params.skip_dedup? Channel.empty() : UMITOOLS_DEDUP_SPLIT.out.dedup_bam
        dedup_bai = params.skip_dedup? Channel.empty() : UMITOOLS_DEDUP_SPLIT.out.dedup_bai
        dedup_log = params.skip_dedup? Channel.empty() : UMITOOLS_DEDUP_SPLIT.out.dedup_log
        dedup_flagstat = params.skip_dedup? Channel.empty() : UMITOOLS_DEDUP_SPLIT.out.dedup_flagstat
        dedup_idxstats = params.skip_dedup? Channel.empty() : UMITOOLS_DEDUP_SPLIT.out.dedup_idxstats

        gene_qc_stats = ch_gene_qc_stats
        transcript_qc_stats = ch_transcript_qc_stats
}
