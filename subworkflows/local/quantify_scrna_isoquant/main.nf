//
// Performs feature quantification for long read single-cell rna data
//

include { ISOQUANT                               } from '../../../modules/local/isoquant'
include { QC_SCRNA as QC_SCRNA_GENE              } from '../../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT        } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_ISOQUANT {
    take:
        in_bam
        _in_bai
        in_flagstat
        in_fasta
        _in_fai
        in_gtf
        _skip_qc
        _skip_seurat

    main:
        ch_versions = channel.empty()

        isoquant_input = in_bam
            .join(in_bai, by: [0])
            .combine(in_fasta.map {_meta, fasta -> [fasta]})
            .combine(in_fai.map {_meta, fai -> [fai]})
            .combine(in_gtf.map {_meta, gtf -> [gtf]})

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
        if (!params.skip_qc && !params.skip_seurat){
            QC_SCRNA_GENE (
                ISOQUANT.out.grouped_gene_mtx_features
                    .join(ISOQUANT.out.grouped_gene_mtx_barcodes, by: [0])
                    .join(ISOQUANT.out.grouped_gene_mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            def new_meta = [ 'id' : meta.id ]
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                in_flagstat,
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
                in_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions                 = ch_versions
        gene_features_file       = ISOQUANT.out.grouped_gene_mtx_features
        gene_barcodes_file       = ISOQUANT.out.grouped_gene_mtx_barcodes
        gene_mtx_file            = ISOQUANT.out.grouped_gene_mtx
        transcript_features_file = ISOQUANT.out.grouped_transcript_mtx_features
        transcript_barcodes_file = ISOQUANT.out.grouped_transcript_mtx_barcodes
        transcript_mtx_file      = ISOQUANT.out.grouped_transcript_mtx
        gene_qc_stats            = ch_gene_qc_stats
        transcript_qc_stats      = ch_transcript_qc_stats
}
