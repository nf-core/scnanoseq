//
// Performs feature quantification for long read single-cell rna data
//

include { ISOQUANT                               } from '../../modules/local/isoquant'
//include { MERGE_MTX as MERGE_MTX_GENE            } from '../../modules/local/merge_mtx'
//include { MERGE_MTX as MERGE_MTX_TRANSCRIPT      } from '../../modules/local/merge_mtx'
//include { SPLIT_GTF                              } from '../../modules/local/split_gtf'
//include { SPLIT_FASTA                            } from '../../modules/local/split_fasta'
//include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SPLIT } from '../../modules/nf-core/samtools/faidx/main'
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

        ISOQUANT (
            in_bam.join(in_bai, by: [0]),
            in_fasta,
            in_fai,
            in_gtf,
            'tag:CB'
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions)

        if (!params.skip_qc && !params.skip_seurat){
            QC_SCRNA_GENE (
                ISOQUANT.out.gene_count_mtx,
                in_flagstat,
                "BASE"
            )
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                ISOQUANT.out.transcript_count_mtx,
                in_flagstat,
                "BASE"
            )
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions = ch_versions
        gene_mtx = ISOQUANT.out.gene_count_mtx
        transcript_mtx = ISOQUANT.out.transcript_count_mtx
        gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
        transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
}
