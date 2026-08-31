//
// Performs alignment-free feature quantification for long read single-cell rna
// data using lr-kallisto. Barcodes and umis come from the flexiplex read names
// rather than from a bam, so this path needs no aligner and no separate
// deduplication step: bustools deduplicates umis while counting.
//
// The reference is prepared once here; the counting itself lives in
// LRKALLISTO_COUNT so it can be run twice over the same index.
//

include { GFFREAD          } from '../../../modules/nf-core/gffread/main'
include { GTF_TO_T2G       } from '../../../modules/local/gtf_to_t2g'
include { KALLISTO_INDEX   } from '../../../modules/local/kallisto/index'
include { LRKALLISTO_COUNT } from '../../../subworkflows/local/lrkallisto_count'

// A second counting pass over the same reads, on XB rather than the assigned
// barcode, so droplets flexiplex could not match to the known list are
// quantified alongside the called cells. Keep the instance above unaliased: the
// selectors in conf/modules.config are keyed on the literal
// '.*:LRKALLISTO_COUNT:...' path.
include { LRKALLISTO_COUNT as LRKALLISTO_COUNT_ALL } from '../../../subworkflows/local/lrkallisto_count'

workflow QUANTIFY_SCRNA_LRKALLISTO {
    take:
        in_fastq          // channel: [ val(meta), path(flexiplex_fastq) ]
        in_known_barcodes // channel: [ val(meta), path(known_barcodes) ]
        in_barcode_counts // channel: [ val(meta), path(barcode_counts) ]
        in_fasta          // channel: [ val(meta), path(genome_fasta) ]
        in_gtf            // channel: [ val(meta), path(gtf) ]
        skip_qc           // bool: Skip qc steps
        skip_seurat       // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        // Barcode and umi widths per assay. Flexiplex trims a fixed-width
        // barcode and umi, which is what lets lr-kallisto read them
        // positionally out of the synthetic barcode read SPLIT_BC_UMI builds.
        def BC_UMI_LENGTHS = [
            '10X_3v3'     : [ 16, 12 ],
            '10X_3v4'     : [ 16, 12 ],
            '10X_5v2'     : [ 16, 10 ],
            '10X_5v3'     : [ 16, 10 ],
            '10X_multiome': [ 16, 12 ],
        ]

        def geometry = BC_UMI_LENGTHS[params.barcode_format]
        if (!geometry) {
            error("lr-kallisto does not support --barcode_format '${params.barcode_format}'. Supported formats: ${BC_UMI_LENGTHS.keySet().join(', ')}")
        }
        def (bc_length, umi_length) = geometry

        // kallisto locates the barcode and umi by position, not by name. File 0
        // is the synthetic barcode+umi read, file 1 is the cdna read.
        def technology = "0,0,${bc_length}:0,${bc_length},${bc_length + umi_length}:1,0,0"

        //
        // MODULE: Extract the transcript sequences to index
        //
        GFFREAD (
            in_gtf,
            in_fasta.map { _meta, fasta -> fasta }
        )

        //
        // MODULE: Build the transcript-to-gene mapping from the same gtf
        //
        GTF_TO_T2G ( in_gtf )
        ch_versions = ch_versions.mix(GTF_TO_T2G.out.versions_gtf_to_t2g)

        //
        // MODULE: Build the lr-kallisto index
        //
        // The genome is only staged when it is going to be used as a d-list.
        // `--kallisto_dlist false` on the command line arrives as the string
        // "false", which is truthy in Groovy, so coerce it rather than testing
        // it directly.
        def use_dlist = params.kallisto_dlist.toString().toBoolean()

        ch_index_input = use_dlist
            ? GFFREAD.out.gffread_fasta.combine( in_fasta.map { _meta, fasta -> [ fasta ] } )
            : GFFREAD.out.gffread_fasta.map { meta, transcripts -> [ meta, transcripts, [] ] }

        KALLISTO_INDEX ( ch_index_input )
        ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions_kallisto_index)

        ch_index = KALLISTO_INDEX.out.index.first()
        ch_t2g = GTF_TO_T2G.out.t2g.first()

        //
        // SUBWORKFLOW: Called cells only. The barcode comes from the read name,
        // which is the one flexiplex matched to the known list, and correcting
        // against that list drops the reads that matched nothing.
        //
        LRKALLISTO_COUNT (
            in_fastq,
            in_known_barcodes,
            in_barcode_counts,
            ch_index,
            ch_t2g,
            technology,
            bc_length,
            umi_length,
            '',
            0,
            true,
            skip_qc,
            skip_seurat
        )
        ch_versions = ch_versions.mix(LRKALLISTO_COUNT.out.versions)

        //
        // SUBWORKFLOW: Every droplet, off the same reads. XB holds the assigned
        // barcode for a called cell and the uncorrected one for a droplet that
        // fell below the knee, so both land in a single barcode space. Only the
        // flexiplex assign module writes XB, which this path already requires:
        // the pipeline refuses --demux_tool_cdna blaze with lr-kallisto. Seurat
        // QC is skipped here, it is a per-cell report and the called cells
        // already have one from the pass above.
        //
        // XB is uncorrected below the knee, so --lrkallisto_all_min_reads puts a
        // floor under the barcode space before the EM ever sees it. Left
        // unbounded a full sample produces millions of single-read barcodes
        // against a few thousand real cells and the EM does not finish.
        if (!params.skip_lrkallisto_all_droplets) {
            LRKALLISTO_COUNT_ALL (
                in_fastq,
                in_known_barcodes,
                in_barcode_counts,
                ch_index,
                ch_t2g,
                technology,
                bc_length,
                umi_length,
                'XB',
                params.lrkallisto_all_min_reads,
                false,
                skip_qc,
                true
            )
            ch_versions = ch_versions.mix(LRKALLISTO_COUNT_ALL.out.versions)
        }

    emit:
        versions                 = ch_versions
        gene_features_file       = LRKALLISTO_COUNT.out.gene_features_file
        gene_barcodes_file       = LRKALLISTO_COUNT.out.gene_barcodes_file
        gene_mtx_file            = LRKALLISTO_COUNT.out.gene_mtx_file
        transcript_features_file = LRKALLISTO_COUNT.out.transcript_features_file
        transcript_barcodes_file = LRKALLISTO_COUNT.out.transcript_barcodes_file
        transcript_mtx_file      = LRKALLISTO_COUNT.out.transcript_mtx_file
        gene_qc_stats            = LRKALLISTO_COUNT.out.gene_qc_stats
        transcript_qc_stats      = LRKALLISTO_COUNT.out.transcript_qc_stats
}
