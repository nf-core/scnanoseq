//
// Performs feature quantification for long read single-cell rna data
//

include { BAMTOOLS_SPLIT                         } from '../../modules/nf-core/bamtools/split/main'
include { ISOQUANT                               } from '../../modules/local/isoquant'
include { MERGE_MTX as MERGE_MTX_GENE            } from '../../modules/local/merge_mtx'
include { MERGE_MTX as MERGE_MTX_TRANSCRIPT      } from '../../modules/local/merge_mtx'
include { QC_SCRNA as QC_SCRNA_GENE              } from '../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT        } from '../../subworkflows/local/qc_scrna'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SPLIT } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SPLIT } from '../../modules/nf-core/samtools/index/main'
include { SPLIT_GTF                              } from '../../modules/local/split_gtf'
include { SPLIT_FASTA                            } from '../../modules/local/split_fasta'

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
        // bam and bai files need to be joined with split fasta, fai and gtf files
        isoquant_input = in_bam
            .join(in_bai, by: [0])
            .combine(in_fasta, by: [0])
            .combine(in_fai, by: [0])
            .combine(in_gtf, by: [0])

        //
        // MODULE: Isoquant
        //
        ISOQUANT (
            isoquant_input,
            'tag:CB'
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions)
    
        if (!params.skip_qc && !params.skip_seurat){
            QC_SCRNA_GENE (
                ISOQUANT.out.grouped_gene_counts,
                in_flagstat,
                "BASE"
            )
            ch_gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                ISOQUANT.out.grouped_transcript_counts,
                in_flagstat,
                "BASE"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions            = ch_versions
        gene_mtx            = Channel.empty()
        transcript_mtx      = Channel.empty()
        gene_qc_stats       = Channel.empty()
        transcript_qc_stats = Channel.empty()
}
