//
// Performs alignment
//

// SUBWORKFLOWS
include { ALIGN_LONGREADS         } from '../../../subworkflows/local/align_longreads'
include { QUANTIFY_SCRNA_ISOQUANT } from '../../../subworkflows/local/quantify_scrna_isoquant'
include { QUANTIFY_SCRNA_OARFISH  } from '../../../subworkflows/local/quantify_scrna_oarfish'
include { DEDUP_UMIS              } from '../../../subworkflows/local/dedup_umis'

// A second pass over the same bam, grouping on XB so that droplets flexiplex could
// not match to the known list are quantified alongside the called cells.
// Keep the instance above unaliased: the QC_SCRNA selectors in conf/modules.config
// are keyed on the literal '.*QUANTIFY_SCRNA_ISOQUANT:QC_SCRNA_...' path.
include { QUANTIFY_SCRNA_ISOQUANT as QUANTIFY_SCRNA_ISOQUANT_ALL } from '../../../subworkflows/local/quantify_scrna_isoquant'

// MODULES
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TAGGED } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TAGGED       } from '../../../modules/nf-core/samtools/index'

include { TAG_BARCODES } from '../../../modules/local/tag_barcodes'


workflow PROCESS_LONGREAD_SCRNA {
    take:
        ch_fasta           // channel: [ val(meta), path(fasta) ]
        ch_fai             // channel: [ val(meta), path(fai) ]
        ch_gtf             // channel: [ val(meta), path(gtf) ]
        ch_fastq           // channel: [ val(meta), path(fastq) ]
        ch_rseqc_bed       // channel: [ val(meta), path(rseqc_bed) ]
        ch_read_bc_info    // channel: [ val(meta), path(read_barcode_info) ]

        val_quant_list      // list: List of quantifiers to use
        val_dedup_tool      // str: Name of deduplication tool to use
        val_genome_aligned  // bool: Whether the bam is aligned to the genome or not
        val_fasta_delimiter // str: Delimiter character used in sequence id in fasta

        val_skip_save_minimap2_index // bool: Skip saving the minimap2 index
        val_skip_qc                  // bool: Skip qc steps
        val_skip_rseqc               // bool: Skip RSeQC
        val_skip_bam_nanocomp        // bool: Skip Nanocomp
        val_skip_seurat              // bool: Skip seurat qc steps
        val_skip_dedup               // bool: Skip deduplication

    main:
        ch_versions = channel.empty()

        //
        // SUBWORKFLOW: Align long Read Data
        //

        ALIGN_LONGREADS(
            ch_fasta,
            ch_fastq,
            ch_rseqc_bed,
            val_skip_save_minimap2_index,
            val_skip_qc,
            val_skip_rseqc,
            val_skip_bam_nanocomp
        )

        //
        // MODULE: Tag Barcodes
        //
        ch_tagged_bam = Channel.empty()
        ch_tagged_bai = Channel.empty()

        if (params.demux_tool_cdna == "flexiplex") {
            //
            // Nothing to do: flexiplex wrote CB/CR/UB/UR into the fastq comment, the
            // assign module derived XB alongside them, and minimap2 -y carried all of
            // it onto the alignments. Reads with no barcode region were dropped during
            // assignment, so every alignment here has a real XB and a full-width UMI --
            // there is nothing umi_tools has to be kept away from.
            //
            ch_tagged_bam = ALIGN_LONGREADS.out.sorted_bam
            ch_tagged_bai = ALIGN_LONGREADS.out.sorted_bai

        } else if (params.demux_tool_cdna == "blaze") {
            TAG_BARCODES (
                ALIGN_LONGREADS.out.sorted_bam
                    .join( ALIGN_LONGREADS.out.sorted_bai, by: 0 )
                    .join( ch_read_bc_info, by: 0)
            )
            ch_versions = ch_versions.mix(TAG_BARCODES.out.versions_tag_barcodes)

            //
            // MODULE: Index Tagged Bam
            //
            SAMTOOLS_INDEX_TAGGED ( TAG_BARCODES.out.tagged_bam )

            ch_tagged_bam = TAG_BARCODES.out.tagged_bam
            ch_tagged_bai = SAMTOOLS_INDEX_TAGGED.out.bai
        }

        //
        // MODULE: Flagstat Tagged Bam
        //
        SAMTOOLS_FLAGSTAT_TAGGED (
            ch_tagged_bam.join(ch_tagged_bai)
        )

        ch_bam = channel.empty()
        ch_bai = channel.empty()
        ch_flagstat = channel.empty()
        ch_idxstats = channel.empty()

        if (!val_skip_dedup) {
            DEDUP_UMIS (
                ch_fasta,
                ch_fai,
                ch_gtf,
                ch_tagged_bam,
                ch_tagged_bai,
                true, // Used to split the bam
                val_genome_aligned,
                val_dedup_tool,
                val_fasta_delimiter
            )

            ch_bam = DEDUP_UMIS.out.dedup_bam
            ch_bai = DEDUP_UMIS.out.dedup_bai
            ch_flagstat = DEDUP_UMIS.out.dedup_flagstat
            ch_versions = DEDUP_UMIS.out.versions
        } else {

            ch_bam = ch_tagged_bam
            ch_bai = ch_tagged_bai
            ch_flagstat = SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat
                .map{
                    meta, flagstat ->
                        def id = meta
                    [id, flagstat]
                }

        }

        //
        // SUBWORKFLOW: Quantify Features
        //

        ch_gene_qc_stats = channel.empty()
        ch_transcript_qc_stats = channel.empty()

        if (val_quant_list.contains("oarfish")) {
            QUANTIFY_SCRNA_OARFISH (
                ch_bam,
                ch_flagstat,
                ch_fasta,
                val_skip_qc,
                val_skip_seurat
            )
            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_OARFISH.out.versions)
            ch_transcript_qc_stats = QUANTIFY_SCRNA_OARFISH.out.transcript_qc_stats
        }

        if (val_quant_list.contains("isoquant")) {
            //
            // Called cells only. Grouping on CB means the reads flexiplex could not
            // match to the known list collect in a single "-" column to be dropped.
            //
            QUANTIFY_SCRNA_ISOQUANT (
                ch_bam,
                ch_bai,
                ch_flagstat,
                ch_fasta,
                ch_fai,
                ch_gtf,
                'tag:CB',
                val_skip_qc,
                val_skip_seurat
            )

            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_ISOQUANT.out.versions)
            ch_gene_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.gene_qc_stats
            ch_transcript_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.transcript_qc_stats

            //
            // Every droplet, off the same alignments. XB holds the corrected barcode
            // for a called cell and the uncorrected one for a droplet that fell below
            // the knee, so both land in a single barcode space. Seurat QC is skipped
            // here: it is a per-cell report, and the called cells already have one
            // from the run above.
            //
            if (params.demux_tool_cdna == "flexiplex") {
                QUANTIFY_SCRNA_ISOQUANT_ALL (
                    ch_bam,
                    ch_bai,
                    ch_flagstat,
                    ch_fasta,
                    ch_fai,
                    ch_gtf,
                    'tag:XB',
                    val_skip_qc,
                    true
                )

                ch_versions = ch_versions.mix(QUANTIFY_SCRNA_ISOQUANT_ALL.out.versions)
            }
        }

    emit:
        // Versions
        versions                 = ch_versions // channel: [ val(meta), path(versions) ]

        // Minimap results + qc's
        minimap_bam              = ALIGN_LONGREADS.out.sorted_bam       // channel: [ val(meta), path(bam) ]
        minimap_bai              = ALIGN_LONGREADS.out.sorted_bai       // channel: [ val(meta), path(bai) ]
        minimap_stats            = ALIGN_LONGREADS.out.stats            // channel: [ val(meta), path(stats) ]
        minimap_flagstat         = ALIGN_LONGREADS.out.flagstat         // channel: [ val(meta), path(flagstat) ]
        minimap_idxstats         = ALIGN_LONGREADS.out.idxstats         // channel: [ val(meta), path(idxstats) ]
        minimap_rseqc_read_dist  = ALIGN_LONGREADS.out.rseqc_read_dist  // channel: [ val(meta), path(rseqc_read_dist) ]
        minimap_nanocomp_bam_txt = ALIGN_LONGREADS.out.nanocomp_bam_txt // channel: [ val(meta), path(nanocomp_bam_txt) ]

        // Barcode tagging results + qc's
        bc_tagged_bam            = ch_tagged_bam                         // channel: [ val(meta), path(bam) ]
        bc_tagged_bai            = ch_tagged_bai                         // channel: [ val(meta), path(bai) ]
        bc_tagged_flagstat       = SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat // channel: [ val(meta), path(flagstat) ]

        // Deduplication results
        dedup_bam                = ch_bam      // channel: [ val(meta), path(bam) ]
        dedup_bai                = ch_bai      // channel: [ val(meta), path(bai) ]
        dedup_flagstat           = ch_flagstat // channel: [ val(meta), path(flagstat) ]
        dedup_idxstats           = ch_idxstats // channel: [ val(meta), path(idxstats) ]

        // Seurat QC Stats
        gene_qc_stats            = ch_gene_qc_stats       // channel: [ val(meta), path(gene_qc_stats) ]
        transcript_qc_stats      = ch_transcript_qc_stats // channel: [ val(meta), path(transcript_qc_stats) ]
}
