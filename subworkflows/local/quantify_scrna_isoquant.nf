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
                    fasta_basename = fasta.toString().split('/')[-1]
                    meta = [ 'chr': fasta_basename.split(/\./)[0] ]
                    [ meta, fasta ]
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
                    gtf_basename = gtf.toString().split('/')[-1]
                    meta = ['chr': gtf_basename.split(/\./)[0]]
                    [ meta, gtf ]
            }
        ch_versions = ch_versions.mix(SPLIT_GTF.out.versions)

        //
        // MODULE: Bamtools split
        //
        BAMTOOLS_SPLIT ( in_bam )
        ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions.first())
        ch_split_bam = BAMTOOLS_SPLIT.out.bam
            .map{
                meta, bam ->
                    [bam]
            }
            .flatten()
            .map{
                bam ->
                    bam_basename = bam.toString().split('/')[-1]
                    split_bam_basename = bam_basename.split(/\./)
                    meta = [
                        'id': split_bam_basename.take(split_bam_basename.size()-1).join("."),
                    ]
                    [ meta, bam ]
            }

        //
        // MODULE: Samtools Index
        //
        SAMTOOLS_INDEX_SPLIT( ch_split_bam )
        ch_split_bai = SAMTOOLS_INDEX_SPLIT.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SPLIT.out.versions.first())

        //
        // MODULE: Isoquant
        //
        ISOQUANT (
            ch_split_bam
                .join(ch_split_bai, by: [0])
                .map{
                    meta, bam, bai ->
                        bam_basename = bam.toString().split('/')[-1]
                        split_bam_basename = bam_basename.split(/\./)
                        chr = [
                            'chr': split_bam_basename[1].replace("REF_","")
                        ]
                        [ chr, meta, bam, bai]
                }
                .combine(ch_split_fasta, by: [0])
                .combine(ch_split_fai, by: [0])
                .combine(ch_split_gtf, by: [0])
                .map{
                    chr, meta, bam, bai, fasta, fai, gtf ->
                        meta.sample_name = meta.id.split(/\./)[0]
                        meta.chr = meta.id.split(/\./)[1]
                        [ meta, bam, bai, fasta, fai, gtf ]
                },
            'tag:CB'
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions)

        //
        // MODULE: Merge Matrix
        //
        MERGE_MTX_GENE (
            ISOQUANT.out.grouped_gene_counts
                .map{
                    meta, mtx ->
                        basename = mtx.toString().split('/')[-1]
                        split_basename = basename.split(/\./)
                        meta = [ 'id': split_basename[0] ]
                    [ meta, mtx ]
                }
                .groupTuple()
        )
        ch_merged_gene_mtx = MERGE_MTX_GENE.out.merged_mtx
        ch_versions = ch_versions.mix(MERGE_MTX_GENE.out.versions)

        MERGE_MTX_TRANSCRIPT (
            ISOQUANT.out.grouped_transcript_counts
                .map{
                    meta, mtx ->
                        basename = mtx.toString().split('/')[-1]
                        split_basename = basename.split(/\./)
                        meta = [ 'id': split_basename[0] ]
                    [ meta, mtx ]
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
