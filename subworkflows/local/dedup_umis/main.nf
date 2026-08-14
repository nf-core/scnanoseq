//
// Rum UMI Deduplication and optionally split the bam for better parallel processing
//

//
// MODULES
//
include { BAMTOOLS_SPLIT                          } from '../../../modules/nf-core/bamtools/split/main'
include { UMITOOLS_DEDUP as UMITOOLS_DEDUP_GENE   } from '../../../modules/nf-core/umitools/dedup/main'
include { UMITOOLS_DEDUP as UMITOOLS_DEDUP_MAIN   } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                          } from '../../../modules/nf-core/samtools/merge/main'
include { SPLIT_BAM                               } from '../../../modules/local/split_bam'
include { GROUP_TRANSCRIPTS                       } from '../../../modules/local/group_transcripts'
include { TAG_GENES                               } from '../../../modules/local/tag_genes'
include { SPLIT_GENE_STATUS                       } from '../../../modules/local/split_gene_status'
include { PICARD_MARKDUPLICATES                   } from '../../../modules/nf-core/picard/markduplicates/main'

//
// SUBWORKFLOWS
//
include { BAM_STATS_SAMTOOLS      } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow DEDUP_UMIS {
    take:
        ch_fasta // channel: [ val(meta), path(fasta) ]
        ch_fai   // channel: [ val(meta), path(fai) ]
        ch_gtf   // channel: [ val(meta), path(gtf) ]
        ch_bam   // channel: [ val(meta), path(bam) ]
        ch_bai   // channel: [ val(meta), path(bai) ]

        val_split_bam       // bool: Split the bam
        val_genome_aligned  // bool: If the bam is aligned to the genome or not
        val_dedup_tool      // str: Name of deduplication tool to use
        val_fasta_delimiter // str: Delimiter character used in the sequence id in fasta

    main:
        ch_versions = channel.empty()

        ch_undedup_bam = channel.empty()
        ch_undedup_bai = channel.empty()

        //
        // Per-gene grouping only makes sense on a genome alignment, and only
        // umi_tools can do it. On a transcriptome alignment the equivalent is
        // --per-contig, which is applied through ext.args and needs no gene tag.
        //
        def val_per_gene = val_genome_aligned && val_dedup_tool == 'umitools' && params.dedup_per_gene

        // Only alignments with an unambiguous gene call are grouped by gene.
        // Reads overlapping two gene bodies, or none, keep the positional
        // result. Keep this in step with --skip-tags-regex in the
        // UMITOOLS_DEDUP_GENE ext.args: the two describe the same split.
        def gene_status_filter = '[GS]=="unique"'

        if (val_per_gene && !params.gtf) {
            error("--dedup_per_gene requires a gtf. Provide --gtf, or disable it with --dedup_per_gene=false.")
        }

        ch_dedup_input_bam = ch_bam

        if (val_per_gene) {
            //
            // MODULE: Tag every alignment with its gene assignment
            //
            // Tagging happens once here rather than once per split so the gtf is
            // only parsed a single time.
            //
            TAG_GENES (
                ch_bam
                    .join(ch_bai)
                    .combine(ch_gtf.map{ meta, gtf -> gtf }.first())
            )
            ch_dedup_input_bam = TAG_GENES.out.gene_tagged_bam
        }

        if (val_split_bam) {
            ch_split_bam = channel.empty()

            if (val_genome_aligned) {
                //
                // MODULE: Bamtools split
                //
                BAMTOOLS_SPLIT ( ch_dedup_input_bam )
                ch_split_bam = BAMTOOLS_SPLIT.out.bam
                    .flatMap{
                        meta, bam ->
                            def bamList = bam instanceof List ? bam : [bam]
                            bamList.collect { b ->
                                def bam_basename = b.toString().split('/')[-1]
                                def split_bam_basename = bam_basename.split(/\./)
                                [ meta + [ 'id': split_bam_basename.take(split_bam_basename.size()-1).join(".") ], b ]
                            }
                    }

            } else {
                //
                // MODULE: Group Transcripts
                //
                GROUP_TRANSCRIPTS (
                    ch_fasta,
                    ch_gtf,
                    val_fasta_delimiter
                )
                ch_versions = ch_versions.mix(GROUP_TRANSCRIPTS.out.versions_group_transcripts)

                //
                // MODULE: Samtools View
                //
                SPLIT_BAM(
                    ch_bam
                        .join(ch_bai)
                        .combine(GROUP_TRANSCRIPTS.out.grouped_transcripts.flatten())
                        .map{
                            meta, bam, bai, region ->
                                def region_basename = region.toString().split('/')[-1]
                                def split_region_basename = region_basename.split(/\./)
                                [['id': meta.id + "." + split_region_basename[0]], bam, bai, region]
                        }
                )
                ch_split_bam = SPLIT_BAM.out.split_bam
                ch_versions = ch_versions.mix(SPLIT_BAM.out.versions_split_bam)
            }

            //
            // SUBWORKFLOW : Sort and index bam
            //
            BAM_SORT_STATS_SAMTOOLS(
                ch_split_bam,
                ch_fasta.first()
            )

            ch_undedup_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_undedup_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

        }
        else {
            ch_undedup_bam = ch_dedup_input_bam
            ch_undedup_bai = ch_bai
        }

        ch_dedup_bam = channel.empty()
        ch_dedup_bai = channel.empty()

        if (val_dedup_tool == 'umitools'){
            if (val_per_gene) {
                //
                // MODULE: Split off the alignments per-gene grouping can act on
                //
                // umi_tools --per-gene drops every read it cannot place on a
                // gene, so the remainder is deduplicated positionally and merged
                // back in below. Without that fallback the run would silently
                // lose around a fifth of its reads.
                //
                SPLIT_GENE_STATUS (
                    ch_undedup_bam.join(ch_undedup_bai, by: [0]),
                    gene_status_filter
                )

                //
                // MODULE: Umitools Dedup, grouping by gene
                //
                UMITOOLS_DEDUP_GENE (
                    SPLIT_GENE_STATUS.out.gene_bam,
                    true )

                //
                // MODULE: Umitools Dedup, positional fallback for the rest
                //
                UMITOOLS_DEDUP_MAIN (
                    SPLIT_GENE_STATUS.out.pos_bam,
                    true )

                ch_dedup_bam = UMITOOLS_DEDUP_GENE.out.bam.mix( UMITOOLS_DEDUP_MAIN.out.bam )

            } else {
                //
                // MODULE: Umitools Dedup
                //
                UMITOOLS_DEDUP_MAIN (
                    ch_undedup_bam.join(ch_undedup_bai, by: [0]),
                    true )
                ch_dedup_bam = UMITOOLS_DEDUP_MAIN.out.bam
            }

        } else {
            //
            // MODULE: Picard Mark Duplicates
            //
            PICARD_MARKDUPLICATES (
                ch_undedup_bam,
                ch_fasta.first(),
                ch_fai.first()
            )
            ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        }

        //
        // MODULE: Samtools Index
        //
        SAMTOOLS_INDEX_DEDUP( ch_dedup_bam )
        ch_dedup_bai = SAMTOOLS_INDEX_DEDUP.out.bai

        if (val_split_bam) {
            //
            // MODULE: Samtools Merge
            //
            SAMTOOLS_MERGE (
                    ch_dedup_bam
                        .map{
                            meta, bam ->
                                def bam_basename = bam.toString().split('/')[-1]
                                def split_bam_basename = bam_basename.split(/\./)
                                def new_meta = meta + [ 'id': split_bam_basename[0] ]
                            [ new_meta, bam ]
                        }
                        .groupTuple()
                        .join(
                            ch_dedup_bai
                                .map{
                                    meta, bai ->
                                        def bai_basename = bai.toString().split('/')[-1]
                                        def split_bai_basename = bai_basename.split(/\./)
                                        def new_meta = meta + [ 'id': split_bai_basename[0] ]
                                    [ new_meta, bai ]
                                }
                                .groupTuple()
                        ),
                ch_fasta
                    .join(ch_fai)
                    .map { meta, fasta_file, fai_file ->
                        [meta, fasta_file, fai_file, "$projectDir/assets/dummy_file.txt"]
                    }
                    .first()

            )
            ch_dedup_bam = SAMTOOLS_MERGE.out.bam

            //
            // MODULE: Samtools Index
            //
            SAMTOOLS_INDEX_MERGED( ch_dedup_bam )
            ch_dedup_bai = SAMTOOLS_INDEX_MERGED.out.bai
        }

        //
        // SUBWORKFLOW: BAM_STATS_SAMTOOLS
        //
        BAM_STATS_SAMTOOLS (
            ch_dedup_bam.join(ch_dedup_bai),
            ch_fasta.first()
        )

    emit:
        versions       = ch_versions                     // channel: [ val(meta), path(versions) ]
        dedup_bam      = ch_dedup_bam                    // channel: [ val(meta), path(bam) ]
        dedup_bai      = ch_dedup_bai                    // channel: [ val(meta), path(bai) ]
        dedup_flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
}
