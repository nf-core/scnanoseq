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

        //
        // MODULE: Split the FASTA
        //
        SPLIT_FASTA( in_fasta )
        ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)
        ch_split_fasta = SPLIT_FASTA.out.split_fasta
            .flatten()
            .map{
                fasta ->
                    def fasta_basename = fasta.toString().split('/')[-1]
                    def new_meta = [ 'chr': fasta_basename.split(/\./)[0] ]
                    [ new_meta, fasta ]
            }

        SAMTOOLS_FAIDX_SPLIT( ch_split_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_split_fai = SAMTOOLS_FAIDX_SPLIT.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SPLIT.out.versions)

        //
        // MODULE: Split the GTF
        //
        SPLIT_GTF( in_gtf )
        ch_split_gtf = SPLIT_GTF.out.split_gtf
            .flatten()
            .map{
                gtf ->
                    def gtf_basename = gtf.toString().split('/')[-1]
                    def new_meta = ['chr': gtf_basename.split(/\./)[0]]
                    [ new_meta, gtf ]
            }
        ch_versions = ch_versions.mix(SPLIT_GTF.out.versions)

        //
        // MODULE: Bamtools split
        //
        BAMTOOLS_SPLIT ( in_bam )
        ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions.first())

        ch_split_bam = BAMTOOLS_SPLIT.out.bam
            .map {
                meta, bam ->
                    return [ bam ]
            }
            .flatten()
            .map { bam ->
                def bam_basename = bam.toString().split('/')[-1]

                def chrom = bam_basename.split(/\./)[1].replace("REF_", "")

                def split_bam_basename = bam_basename.split(/\./)
                def new_id = split_bam_basename.take(split_bam_basename.size()-2).join(".")

                def new_meta = ['id': new_id, 'chr': chrom]
                return [ new_meta, bam ]
            }

        //
        // MODULE: Samtools Index
        //
        SAMTOOLS_INDEX_SPLIT( ch_split_bam )
        ch_split_bai = SAMTOOLS_INDEX_SPLIT.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SPLIT.out.versions.first())

        // Prepare isoquant input channel
        // bam and bai files need to be joined with split fasta, fai and gtf files
        isoquant_input = ch_split_bam
            .join(ch_split_bai, by: [0])
            .map { meta, bam, bai ->
                def chrom = [ 'chr': meta.chr ]
                [ chrom, meta, bam, bai ]
            }
            .combine(ch_split_fasta, by: [0])
            .combine(ch_split_fai, by: [0])
            .combine(ch_split_gtf, by: [0])
            .map { chrom, meta, bam, bai, fasta, fai, gtf ->
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

        //
        // MODULE: Merge Matrix
        //
        MERGE_MTX_GENE (
            ISOQUANT.out.grouped_gene_counts
                .groupTuple()
        )
        ch_merged_gene_mtx = MERGE_MTX_GENE.out.merged_mtx
        ch_versions = ch_versions.mix(MERGE_MTX_GENE.out.versions)

        MERGE_MTX_TRANSCRIPT (
            ISOQUANT.out.grouped_transcript_counts
                .map{
                    meta, mtx ->
                        def basename = mtx.toString().split('/')[-1]
                        def split_basename = basename.split(/\./)
                        def new_meta = [ 'id': split_basename[0] ]
                    [ new_meta, mtx ]
                }
                .groupTuple()
        )
        ch_merged_transcript_mtx = MERGE_MTX_TRANSCRIPT.out.merged_mtx
        ch_versions = ch_versions.mix(MERGE_MTX_TRANSCRIPT.out.versions)

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
