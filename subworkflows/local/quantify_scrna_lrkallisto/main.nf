//
// Performs alignment-free feature quantification for long read single-cell rna
// data using lr-kallisto. Barcodes and umis come from the flexiplex read names
// rather than from a bam, so this path needs no aligner and no separate
// deduplication step: bustools deduplicates umis while counting.
//

include { GFFREAD                         } from '../../../modules/nf-core/gffread/main'
include { GTF_TO_T2G                      } from '../../../modules/local/gtf_to_t2g'
include { KALLISTO_INDEX                  } from '../../../modules/local/kallisto/index'
include { SPLIT_BC_UMI                    } from '../../../modules/local/split_bc_umi'
include { KALLISTO_BUS                    } from '../../../modules/local/kallisto/bus'
include { BUSTOOLS_TCC                    } from '../../../modules/local/bustools/tcc'
include { KALLISTO_QUANTTCC               } from '../../../modules/local/kallisto/quanttcc'
include { QC_SCRNA as QC_SCRNA_GENE       } from '../../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_LRKALLISTO {
    take:
        in_fastq          // channel: [ val(meta), path(flexiplex_fastq) ]
        in_known_barcodes // channel: [ val(meta), path(known_barcodes) ]
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
            '10X_5v3'     : [ 16, 12 ],
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
        // MODULE: Move the barcode and umi out of the read name into their own fastq
        //
        SPLIT_BC_UMI ( in_fastq, bc_length, umi_length )
        ch_versions = ch_versions.mix(SPLIT_BC_UMI.out.versions_split_bc_umi)

        //
        // MODULE: Pseudoalign the reads
        //
        KALLISTO_BUS ( SPLIT_BC_UMI.out.reads, ch_index, technology )
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
        BUSTOOLS_TCC ( ch_tcc_input, ch_t2g )
        ch_versions = ch_versions.mix(BUSTOOLS_TCC.out.versions_bustools_tcc)

        //
        // MODULE: Quantify transcript and gene abundances with the long read EM
        //
        KALLISTO_QUANTTCC (
            BUSTOOLS_TCC.out.tcc_mtx
                .join( BUSTOOLS_TCC.out.tcc_ec, by: [0] )
                .join( BUSTOOLS_TCC.out.barcodes, by: [0] ),
            ch_index,
            ch_t2g
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
