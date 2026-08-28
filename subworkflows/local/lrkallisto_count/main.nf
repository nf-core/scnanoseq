//
// One lr-kallisto counting pass: barcodes and umis out of the flexiplex read,
// pseudoalignment, umi deduplication and the long read EM.
//
// The pass runs twice off the same reads and the same index, once over the
// barcodes flexiplex matched to the known list and once over every droplet.
// Which barcode it counts on is decided by val_barcode_tag and val_correct.
//

include { SPLIT_BC_UMI                    } from '../../../modules/local/split_bc_umi'
include { KALLISTO_BUS                    } from '../../../modules/local/kallisto/bus'
include { BUSTOOLS_TCC                    } from '../../../modules/local/bustools/tcc'
include { KALLISTO_QUANTTCC               } from '../../../modules/local/kallisto/quanttcc'
include { QC_SCRNA as QC_SCRNA_GENE       } from '../../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT } from '../../../subworkflows/local/qc_scrna'

workflow LRKALLISTO_COUNT {
    take:
        in_fastq          // channel: [ val(meta), path(flexiplex_fastq) ]
        in_known_barcodes // channel: [ val(meta), path(known_barcodes) ]
        in_index          // channel: [ val(meta), path(kallisto_index) ]
        in_t2g            // channel: [ val(meta), path(t2g) ]
        val_technology    // str: The kallisto -x geometry string
        val_bc_length     // int: Barcode width
        val_umi_length    // int: Umi width
        val_barcode_tag   // str: Fastq comment tag to take the barcode from, '' for the read name
        val_correct       // bool: Correct barcodes against the known-barcode list
        skip_qc           // bool: Skip qc steps
        skip_seurat       // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        //
        // MODULE: Move the barcode and umi out of the read name into their own fastq
        //
        SPLIT_BC_UMI ( in_fastq, val_bc_length, val_umi_length, val_barcode_tag )
        ch_versions = ch_versions.mix(SPLIT_BC_UMI.out.versions_split_bc_umi)

        //
        // MODULE: Pseudoalign the reads
        //
        KALLISTO_BUS ( SPLIT_BC_UMI.out.reads, in_index, val_technology )
        ch_versions = ch_versions.mix(KALLISTO_BUS.out.versions_kallisto_bus)

        // The demultiplexing subworkflow tags the reads and the known-barcode
        // list with different metadata, so match on the sample id rather than
        // joining on the whole map.
        ch_tcc_input = KALLISTO_BUS.out.bus
            .join( KALLISTO_BUS.out.ecmap, by: [0] )
            .join( KALLISTO_BUS.out.txnames, by: [0] )
            .combine( in_known_barcodes )
            .filter { meta, _bus, _ec, _tx, meta2, _barcodes -> meta.id == meta2.id }
            .map { meta, bus, ec, tx, _meta2, barcodes -> [ meta, bus, ec, tx, barcodes ] }

        //
        // MODULE: Correct barcodes, deduplicate umis and count equivalence classes
        //
        BUSTOOLS_TCC ( ch_tcc_input, in_t2g, val_correct )
        ch_versions = ch_versions.mix(BUSTOOLS_TCC.out.versions_bustools_tcc)

        //
        // MODULE: Quantify transcript and gene abundances with the long read EM
        //
        KALLISTO_QUANTTCC (
            BUSTOOLS_TCC.out.tcc_mtx
                .join( BUSTOOLS_TCC.out.tcc_ec, by: [0] )
                .join( BUSTOOLS_TCC.out.barcodes, by: [0] ),
            in_index,
            in_t2g
        )
        ch_versions = ch_versions.mix(KALLISTO_QUANTTCC.out.versions_kallisto_quanttcc)

        ch_gene_qc_stats = channel.empty()
        ch_transcript_qc_stats = channel.empty()
        if (!skip_qc && !skip_seurat) {

            ch_flagstat = SPLIT_BC_UMI.out.flagstat
                .map { meta, flagstat -> [ [ 'id': meta.id, 'type': meta.type ], flagstat ] }

            QC_SCRNA_GENE (
                KALLISTO_QUANTTCC.out.gene_features
                    .join( KALLISTO_QUANTTCC.out.gene_barcodes, by: [0] )
                    .join( KALLISTO_QUANTTCC.out.gene_mtx, by: [0] )
                    .map { meta, features, barcodes, mtx ->
                        [ [ 'id': meta.id, 'type': meta.type ], [ features, barcodes, mtx ] ]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                KALLISTO_QUANTTCC.out.transcript_features
                    .join( KALLISTO_QUANTTCC.out.transcript_barcodes, by: [0] )
                    .join( KALLISTO_QUANTTCC.out.transcript_mtx, by: [0] )
                    .map { meta, features, barcodes, mtx ->
                        [ [ 'id': meta.id, 'type': meta.type ], [ features, barcodes, mtx ] ]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions                 = ch_versions
        gene_features_file       = KALLISTO_QUANTTCC.out.gene_features
        gene_barcodes_file       = KALLISTO_QUANTTCC.out.gene_barcodes
        gene_mtx_file            = KALLISTO_QUANTTCC.out.gene_mtx
        transcript_features_file = KALLISTO_QUANTTCC.out.transcript_features
        transcript_barcodes_file = KALLISTO_QUANTTCC.out.transcript_barcodes
        transcript_mtx_file      = KALLISTO_QUANTTCC.out.transcript_mtx
        gene_qc_stats            = ch_gene_qc_stats
        transcript_qc_stats      = ch_transcript_qc_stats
}
