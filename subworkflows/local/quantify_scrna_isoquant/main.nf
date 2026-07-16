//
// Performs feature quantification for long read single-cell rna data
//

include { ISOQUANT                               } from '../../../modules/local/isoquant'
include { QC_SCRNA as QC_SCRNA_GENE              } from '../../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT        } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_ISOQUANT {
    take:
        ch_bam      // channel: [ val(meta), path(bam) ]
        ch_bai      // channel: [ val(meta), path(bai) ]
        ch_flagstat // channel: [ val(meta), path(flagstat) ]
        ch_fasta    // channel: [ val(meta), path(fasta) ]
        ch_fai      // channel: [ val(meta), path(fai) ]
        ch_gtf      // channel: [ val(meta), path(gtf) ]

        val_skip_qc     // bool: Skip qc steps
        val_skip_seurat // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        isoquant_input = ch_bam
            .join(ch_bai, by: [0])
            .combine(ch_fasta.map {_meta, fasta -> [fasta]})
            .combine(ch_fai.map {_meta, fai -> [fai]})
            .combine(ch_gtf.map {_meta, gtf -> [gtf]})

        //
        // MODULE: Isoquant
        //
        ISOQUANT (
            isoquant_input,
            'tag:CB'
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions_isoquant)

        ch_gene_qc_stats = channel.empty()
        ch_transcript_qc_stats = channel.empty()
        if (!val_skip_qc && !val_skip_seurat){
            QC_SCRNA_GENE (
                ISOQUANT.out.grouped_gene_mtx_features
                    .join(ISOQUANT.out.grouped_gene_mtx_barcodes, by: [0])
                    .join(ISOQUANT.out.grouped_gene_mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            def new_meta = [ 'id' : meta.id ]
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                ISOQUANT.out.grouped_transcript_mtx_features
                    .join(ISOQUANT.out.grouped_transcript_mtx_barcodes, by: [0])
                    .join(ISOQUANT.out.grouped_transcript_mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            def new_meta = [ 'id' : meta.id ]
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions                 = ch_versions                                  // channel: [ val(meta), path(versions) ]
        gene_features_file       = ISOQUANT.out.grouped_gene_mtx_features       // channel: [ val(meta), path(gene_features) ]
        gene_barcodes_file       = ISOQUANT.out.grouped_gene_mtx_barcodes       // channel: [ val(meta), path(gene_barcodes) ]
        gene_mtx_file            = ISOQUANT.out.grouped_gene_mtx                // channel: [ val(meta), path(gene_mtx) ]
        transcript_features_file = ISOQUANT.out.grouped_transcript_mtx_features // channel: [ val(meta), path(transcript_features) ]
        transcript_barcodes_file = ISOQUANT.out.grouped_transcript_mtx_barcodes // channel: [ val(meta), path(transcript_barcodes) ]
        transcript_mtx_file      = ISOQUANT.out.grouped_transcript_mtx          // channel: [ val(meta), path(transcript_mtx) ]
        gene_qc_stats            = ch_gene_qc_stats                             // channel: [ val(meta), path(gene_qc_stats) ]
        transcript_qc_stats      = ch_transcript_qc_stats                       // channel: [ val(meta), path(transcript_qc_stats) ]
}
