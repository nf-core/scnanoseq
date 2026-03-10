//
// Performs feature quantification for long read single-cell rna data
//

include { ISOQUANT                               } from '../../modules/local/isoquant'
include { QC_SCRNA as QC_SCRNA_GENE              } from '../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT        } from '../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_ISOQUANT {
    take:
        in_bam
        in_bai
        in_flagstat
        in_fasta
        in_fai
        in_gtf
        skip_qc
        skip_seurat

    main:
        ch_versions = Channel.empty()

        // Prepare isoquant input channel
        // Run IsoQuant once per sample on full BAM + full reference files.
        ch_fasta = in_fasta.map { meta, fasta -> fasta }
        ch_fai   = in_fai.map   { meta, fai   -> fai }
        ch_gtf   = in_gtf.map   { meta, gtf   -> gtf }

        isoquant_input = in_bam
            .join(in_bai, by: [0])
            .combine(ch_fasta)
            .combine(ch_fai)
            .combine(ch_gtf)
            .map { meta, bam, bai, fasta, fai, gtf ->
                [ meta, bam, bai, fasta, fai, gtf ]
            }

        //
        // MODULE: Isoquant
        //

        ISOQUANT (
            isoquant_input,
            'tag:CB'
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions)

        // Use IsoQuant grouped outputs directly (single run per sample, no split/merge).
        ch_merged_gene_mtx = ISOQUANT.out.grouped_gene_counts
            .map {
                meta, gene_mtx ->
                    def new_meta = [ 'id': meta.id ]
                    return [ new_meta, gene_mtx ]
            }

        ch_merged_transcript_mtx = ISOQUANT.out.grouped_transcript_counts
            .map {
                meta, transcript_mtx ->
                    def new_meta = [ 'id': meta.id ]
                    return [ new_meta, transcript_mtx ]
            }

        ch_gene_qc_stats = Channel.empty()
        ch_transcript_qc_stats = Channel.empty()

        if (!params.skip_qc && !params.skip_seurat){
            QC_SCRNA_GENE (
                ch_merged_gene_mtx,
                in_flagstat,
                "BASE"
            )
            ch_gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                ch_merged_transcript_mtx,
                in_flagstat,
                "BASE"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions            = ch_versions
        gene_mtx            = ch_merged_gene_mtx
        transcript_mtx      = ch_merged_transcript_mtx
        gene_qc_stats       = ch_gene_qc_stats
        transcript_qc_stats = ch_transcript_qc_stats
}
