//
// Performs alignment
//

// SUBWORKFLOWS
include { ALIGN_LONGREADS         } from '../../../subworkflows/local/align_longreads'
include { QUANTIFY_SCRNA_ISOQUANT } from '../../../subworkflows/local/quantify_scrna_isoquant'
include { QUANTIFY_SCRNA_OARFISH  } from '../../../subworkflows/local/quantify_scrna_oarfish'
include { DEDUP_UMIS              } from '../../../subworkflows/local/dedup_umis'

// A second pass over the all-reads bam, grouping on XB so that droplets flexiplex
// could not match to the known list are quantified alongside the called cells.
// Keep the instance above unaliased: the QC_SCRNA selectors in conf/modules.config
// are keyed on the literal '.*QUANTIFY_SCRNA_ISOQUANT:QC_SCRNA_...' path.
include { QUANTIFY_SCRNA_ISOQUANT as QUANTIFY_SCRNA_ISOQUANT_ALL } from '../../../subworkflows/local/quantify_scrna_isoquant'

// MODULES
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TAGGED } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DEDUP  } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TAGGED       } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP        } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_DEDUP        } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_ALL_READS    } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALL_READS    } from '../../../modules/nf-core/samtools/index'

include { TAG_BARCODES         } from '../../../modules/local/tag_barcodes'
include { SPLIT_BARCODE_STATUS } from '../../../modules/local/split_barcode_status'


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

        // Alignments that skip deduplication and are merged back in afterwards.
        ch_undedupable_bam = Channel.empty()

        // Reads flexiplex could place in no droplet at all. Keep this in step with
        // the XB rule in modules/local/flexiplex/assign: XB falls back to "-" only
        // when neither a corrected nor an uncorrected barcode was found.
        def barcode_status_filter = '[XB]!="-"'

        if (params.demux_tool_cdna == "flexiplex") {
            //
            // MODULE: Split off the alignments UMI deduplication can act on
            //
            // No tagging step is needed: flexiplex wrote CB/CR/UB/UR into the fastq
            // comment, the pipeline derived XB alongside them, and minimap2 -y carried
            // all of it onto the alignments. What is needed is to keep the reads with
            // no barcode away from umi_tools, which has no UMI to group them by.
            //
            SPLIT_BARCODE_STATUS (
                ALIGN_LONGREADS.out.sorted_bam.join( ALIGN_LONGREADS.out.sorted_bai, by: 0 ),
                barcode_status_filter
            )

            ch_tagged_bam      = SPLIT_BARCODE_STATUS.out.barcoded_bam.map { meta, bam, _bai -> [ meta, bam ] }
            ch_tagged_bai      = SPLIT_BARCODE_STATUS.out.barcoded_bam.map { meta, _bam, bai -> [ meta, bai ] }
            ch_undedupable_bam = SPLIT_BARCODE_STATUS.out.nobarcode_bam.map { meta, bam, _bai -> [ meta, bam ] }

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
        // MODULE: Merge the reads that bypassed deduplication back in
        //
        // The published bam is meant to hold every read, so the unbarcoded alignments
        // rejoin the deduplicated ones here. They are merged rather than deduplicated
        // because they carry no UMI to group on. ch_bam is deliberately left pointing
        // at the deduplicated, barcoded alignments so that everything keyed on a called
        // cell -- oarfish and the filtered isoquant run -- sees what it sees today.
        //
        ch_all_reads_bam      = ch_bam
        ch_all_reads_bai      = ch_bai
        ch_all_reads_flagstat = ch_flagstat

        if (params.demux_tool_cdna == "flexiplex") {
            SAMTOOLS_MERGE_ALL_READS (
                ch_bam
                    .mix( ch_undedupable_bam )
                    .map { meta, bam ->
                        [ meta.subMap('id', 'single_end', 'cell_count', 'type'), bam ] }
                    .groupTuple()
                    .map { meta, bams -> [ meta, bams, [] ] },
                ch_fasta
                    .join(ch_fai)
                    .map { meta, fasta_file, fai_file ->
                        [ meta, fasta_file, fai_file, [] ] }
                    .first()
            )

            SAMTOOLS_INDEX_ALL_READS ( SAMTOOLS_MERGE_ALL_READS.out.bam )

            SAMTOOLS_FLAGSTAT_DEDUP (
                SAMTOOLS_MERGE_ALL_READS.out.bam.join( SAMTOOLS_INDEX_ALL_READS.out.bai )
            )

            ch_all_reads_bam      = SAMTOOLS_MERGE_ALL_READS.out.bam
            ch_all_reads_bai      = SAMTOOLS_INDEX_ALL_READS.out.bai
            ch_all_reads_flagstat = SAMTOOLS_FLAGSTAT_DEDUP.out.flagstat
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
            // Every droplet. XB holds the corrected barcode for a called cell and the
            // uncorrected one for a droplet that fell below the knee, so both land in a
            // single barcode space. Seurat QC is skipped here: it is a per-cell report,
            // and the called cells already have one from the run above.
            //
            if (params.demux_tool_cdna == "flexiplex") {
                QUANTIFY_SCRNA_ISOQUANT_ALL (
                    ch_all_reads_bam,
                    ch_all_reads_bai,
                    ch_all_reads_flagstat,
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
