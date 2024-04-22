/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER PRESETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This is if the user passes in direct regex
def cell_barcode_pattern = ""

// This is for if the user wants to do more human readable regex
def cell_barcode_lengths = ""
def umi_lengths = ""
def blaze_whitelist = ""

// TODO: Move this to a config file
if (params.barcode_format = "cellranger_3_prime") {
    blaze_whitelist = file("$baseDir/assets/whitelist/3M-february-2018.zip")
    bc_length = "16"
    umi_length = "12"
}

if (params.whitelist) {
    blaze_whitelist = params.whitelist
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { NANOFILT                                                 } from "../modules/local/nanofilt"
include { NANOCOMP as NANOCOMP_FASTQ                               } from "../modules/local/nanocomp"
include { NANOCOMP as NANOCOMP_BAM                                 } from "../modules/local/nanocomp"
include { SPLIT_FILE                                               } from "../modules/local/split_file"
include { PIGZ as ZIP_TRIM                                         } from "../modules/local/pigz"
include { BLAZE                                                    } from "../modules/local/blaze"
include { PREEXTRACT_FASTQ                                         } from "../modules/local/preextract_fastq.nf"
include { PAFTOOLS                                                 } from "../modules/local/paftools"
include { MINIMAP2_INDEX                                           } from "../modules/local/minimap2_index"
include { MINIMAP2_ALIGN                                           } from "../modules/local/minimap2_align"
include { TAG_BARCODES                                             } from "../modules/local/tag_barcodes"
include { CORRECT_BARCODES                                         } from "../modules/local/correct_barcodes"
include { ISOQUANT                                                 } from "../modules/local/isoquant"
include { SEURAT as SEURAT_GENE                                    } from "../modules/local/seurat"
include { SEURAT as SEURAT_TRANSCRIPT                              } from "../modules/local/seurat"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_GENE        } from "../modules/local/combine_seurat_stats"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_TRANSCRIPT  } from "../modules/local/combine_seurat_stats"
include { UCSC_GTFTOGENEPRED                                       } from "../modules/local/ucsc_gtftogenepred"
include { UCSC_GENEPREDTOBED                                       } from "../modules/local/ucsc_genepredtobed"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_REFERENCE_FILES } from "../subworkflows/local/prepare_reference_files"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP                                } from "../modules/nf-core/gunzip/main"
include { MULTIQC as MULTIQC_RAWQC              } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_FINALQC            } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { UMITOOLS_DEDUP                        } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BAM    } from "../modules/nf-core/samtools/view/main"
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER } from "../modules/nf-core/samtools/view/main"
include { CAT_CAT                               } from "../modules/nf-core/cat/cat/main"
include { CAT_FASTQ                             } from '../modules/nf-core/cat/fastq/main'
include { RSEQC_READDISTRIBUTION                } from '../modules/nf-core/rseqc/readdistribution/main'
include { BAMTOOLS_SPLIT                        } from '../modules/nf-core/bamtools/split/main'
include { SAMTOOLS_MERGE                        } from '../modules/nf-core/samtools/merge/main'
include { paramsSummaryMap } from 'plugin/nf-validation'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM         } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM        } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT     } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_MINIMAP  } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FILTERED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_TAGGED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_CORRECTED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_SPLIT } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_DEDUP } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scnanoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCNANOSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_samplesheet
        .branch{
            meta, fastq ->
                single: fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastqs }

    //
    // MODULE: Combine fastqs from the same sample
    //
    CAT_FASTQ ( ch_fastqs.multiple )
        .reads
        .mix ( ch_fastqs.single )
        .set { ch_cat_fastq }

    ch_versions = ch_versions.mix (CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot, ToulligQC and FastQC - pre-trim QC
    //

    ch_fastqc_multiqc_pretrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_PRE_TRIM ( ch_cat_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = Channel.empty()
    ch_nanocomp_fastq_txt = Channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ ( ch_cat_fastq.collect{it[1]},
                        ch_dummy_file )
        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.txt

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

    }

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES ( "",
                            "",
                            params.fasta,
                            params.gtf)

    fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    fai = PREPARE_REFERENCE_FILES.out.prepped_fai
    gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf


    ch_versions = ch_versions.mix( PREPARE_REFERENCE_FILES.out.versions )

    //
    // MODULE: Generate junction file - paftools
    //

    PAFTOOLS ( gtf.map { meta, gtf -> [gtf]} )
    ch_bed = PAFTOOLS.out.bed
    ch_versions = ch_versions.mix(PAFTOOLS.out.versions)

    //
    // MODULE: Generate bed file from input gtf for rseqc
    //

    //TODO: This uses params.gtf instead of gtf in PAFTOOLS
    // come back to this once intron work is finished (likely input will be fine)
    ch_pred = Channel.empty()
    ch_rseqc_bed = Channel.empty()
    if (!params.skip_qc && !params.skip_rseqc) {
        UCSC_GTFTOGENEPRED( params.gtf )
        ch_pred = UCSC_GTFTOGENEPRED.out.genepred
        ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

        UCSC_GENEPREDTOBED ( ch_pred )
        ch_rseqc_bed = UCSC_GENEPREDTOBED.out.bed
        ch_versions = ch_versions.mix(UCSC_GENEPREDTOBED.out.versions)
    }

    //
    // MODULE: Trim and filter reads
    //
    ch_zipped_reads = Channel.empty()
    ch_fastqc_multiqc_postrim = Channel.empty()

    if (!params.skip_trimming){
        // The trimmers require unzipped fastqs, so we'll need to unzip them
        // to accomodate this requirement

        //
        // MODULE: Unzip fastq
        //
        GUNZIP( ch_cat_fastq )
        ch_unzipped_fastqs = GUNZIP.out.gunzip
        ch_versions = ch_versions.mix( GUNZIP.out.versions )

        //
        // MODULE: Split fastq
        //
        ch_fastqs = ch_unzipped_fastqs

        if (params.split_amount > 0) {
            SPLIT_FILE( ch_unzipped_fastqs, '.fastq', params.split_amount )

            // Temporarily change the meta object so that the id is present on the
            // fastq to prevent duplicated names
            SPLIT_FILE.out.split_files
                .transpose()
                .set { ch_fastqs }

            ch_versions = ch_versions.mix(SPLIT_FILE.out.versions)
        }

        ch_trimmed_reads = ch_fastqs
        if (!params.skip_trimming) {

            NANOFILT ( ch_fastqs )
            ch_trimmed_reads = NANOFILT.out.reads
            ch_versions = ch_versions.mix(NANOFILT.out.versions)
        }

        // If the fastqs were split, combine them together
        ch_trimmed_reads_combined = ch_trimmed_reads
        if (params.split_amount > 0){
            CAT_CAT(ch_trimmed_reads.groupTuple())
            ch_trimmed_reads_combined = CAT_CAT.out.file_out
        }

        //
        // MODULE: Zip the reads
        //
        ZIP_TRIM (ch_trimmed_reads_combined, "filtered" )
        ch_zipped_reads = ZIP_TRIM.out.archive
        ch_versions = ch_versions.mix(ZIP_TRIM.out.versions)

        //
        // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
        //
        if (!params.skip_qc){

            //
            // MODULE: Run qc on the post trimmed reads
            //
            FASTQC_NANOPLOT_POST_TRIM ( ch_zipped_reads, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.toulligqc_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.fastqc_version.first().ifEmpty(null))
        }
    } else {
        ch_zipped_reads = ch_cat_fastq
    }

    //
    // MODULE: Generate whitelist
    //

    BLAZE ( ch_zipped_reads, blaze_whitelist)

    ch_putative_bc = BLAZE.out.putative_bc
    ch_gt_whitelist = BLAZE.out.whitelist
    ch_whitelist_bc_count = BLAZE.out.bc_count
    ch_versions = ch_versions.mix(BLAZE.out.versions)

    //
    // MODULE: Extract barcodes
    //

    PREEXTRACT_FASTQ( ch_zipped_reads.join(ch_putative_bc), params.barcode_format)
    ch_zipped_r1_reads = PREEXTRACT_FASTQ.out.r1_reads
    ch_zipped_r2_reads = PREEXTRACT_FASTQ.out.r2_reads

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_EXTRACT ( ch_zipped_r2_reads, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_version.first().ifEmpty(null))
    }

    //
    // MINIMAP2_INDEX
    //

    if (!params.skip_save_minimap2_index) {
        MINIMAP2_INDEX ( fasta.map { meta, fasta -> [fasta]},  ch_bed)
        ch_minimap_index = MINIMAP2_INDEX.out.index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MINIMAP2_ALIGN
    //

    if (!params.skip_save_minimap2_index) {
        ch_reference = ch_minimap_index.toList()
    } else {
        ch_reference = Channel.fromPath(fasta, checkIfExists: true).toList()
    }
    MINIMAP2_ALIGN ( ch_zipped_r2_reads, ch_bed, ch_reference )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.sam
        .combine( ch_dummy_file )
        .set { ch_minimap_sam }

    //
    // MODULE: Samtools view
    //
    SAMTOOLS_VIEW_BAM ( ch_minimap_sam, [[],[]], [] )

    ch_minimap_bam = SAMTOOLS_VIEW_BAM.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_BAM.out.versions)

    // acquire only mapped reads from bam for downstream processing
    // NOTE: some QCs steps are performed on the full BAM
    ch_minimap_bam
        .combine( ch_dummy_file )
        .set { ch_minimap_bam_filter }

    SAMTOOLS_VIEW_FILTER ( ch_minimap_bam_filter, [[],[]], [] )
    ch_minimap_mapped_only_bam = SAMTOOLS_VIEW_FILTER.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILTER.out.versions)

    //
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
    BAM_SORT_STATS_SAMTOOLS_MINIMAP ( ch_minimap_bam,
                                        fasta )
    ch_minimap_sorted_bam = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bam
    ch_minimap_sorted_bai = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bai

    // these stats go for multiqc
    ch_minimap_sorted_stats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.stats
    ch_minimap_sorted_flagstat = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.flagstat
    ch_minimap_sorted_idxstats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.idxstats
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.versions)

    BAM_SORT_STATS_SAMTOOLS_FILTERED ( ch_minimap_mapped_only_bam,
                                        fasta )
    ch_minimap_filtered_sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam
    ch_minimap_filtered_sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_FILTERED.out.versions)

    //
    // MODULE: RSeQC read distribution for BAM files (unfiltered for QC purposes)
    //
    ch_rseqc_read_dist = Channel.empty()
    if (!params.skip_qc && !params.skip_rseqc) {
        RSEQC_READDISTRIBUTION ( ch_minimap_sorted_bam, ch_rseqc_bed )
        ch_rseqc_read_dist = RSEQC_READDISTRIBUTION.out.txt
        ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions)
    }

    //
    // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
    //
    ch_nanocomp_bam_html = Channel.empty()
    ch_nanocomp_bam_txt = Channel.empty()

    if (!params.skip_qc && !params.skip_bam_nanocomp) {

        NANOCOMP_BAM ( ch_minimap_sorted_bam.collect{it[1]},
                        ch_minimap_sorted_bai.collect{it[1]})

        ch_nanocomp_bam_html = NANOCOMP_BAM.out.html
        ch_nanocomp_bam_txt = NANOCOMP_BAM.out.txt
        ch_versions = ch_versions.mix( NANOCOMP_BAM.out.versions )
    }


    //
    // MODULE: Tag Barcodes
    //

    TAG_BARCODES (
        ch_minimap_filtered_sorted_bam
            .join( ch_minimap_filtered_sorted_bai, by: 0)
            .join( ch_zipped_r1_reads, by: 0 ),
        bc_length,
        umi_length
    )

    ch_tagged_bam = TAG_BARCODES.out.tagged_bam
    ch_versions = ch_versions.mix(TAG_BARCODES.out.versions)

    BAM_SORT_STATS_SAMTOOLS_TAGGED ( ch_tagged_bam,
                                        fasta )

    ch_tagged_sorted_bam = BAM_SORT_STATS_SAMTOOLS_TAGGED.out.bam
    ch_tagged_sorted_bai = BAM_SORT_STATS_SAMTOOLS_TAGGED.out.bai

    //
    // MODULE: Correct Barcodes
    //

    CORRECT_BARCODES (
        ch_tagged_sorted_bam
            .join ( ch_tagged_sorted_bai, by: 0)
            .join ( ch_gt_whitelist, by: 0)
            .join ( ch_whitelist_bc_count, by: 0 )
    )

    ch_corrected_bam = CORRECT_BARCODES.out.corrected_bam
    ch_versions = ch_versions.mix(CORRECT_BARCODES.out.versions)

    //
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
    BAM_SORT_STATS_SAMTOOLS_CORRECTED ( ch_corrected_bam,
                                        fasta )

    ch_corrected_sorted_bam = BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.bam
    ch_corrected_sorted_bai = BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.bai

    // these stats go for multiqc
    ch_corrected_sorted_stats =     BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.stats
    ch_corrected_sorted_flagstat =  BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.flagstat
    ch_corrected_sorted_idxstats =  BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.idxstats
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_CORRECTED.out.versions)

    // TODO: Rename the dedup_bam channel to be more descriptive
    ch_dedup_sorted_bam = ch_corrected_sorted_bam
    ch_dedup_sorted_bam_bai = ch_corrected_sorted_bai
    ch_dedup_sorted_flagstat = ch_corrected_sorted_flagstat
    ch_dedup_log = Channel.empty()

    if (!params.skip_dedup) {

        //
        // MODULE: Bamtools Split
        //
        BAMTOOLS_SPLIT ( ch_corrected_sorted_bam )
        ch_split_bams = BAMTOOLS_SPLIT.out.bam

        ch_split_corrected_bam = ch_split_bams
                                    .map{
                                        meta, bam ->
                                            [bam]
                                    }
                                    .flatten()
                                    .map{
                                        bam ->
                                            bam_basename = bam.toString().split('/')[-1]
                                            split_bam_basename = bam_basename.split(/\./)
                                            meta = [ 'id': split_bam_basename.take(split_bam_basename.size()-1).join(".") ]
                                            [ meta, bam ]
                                    }

        //
        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        BAM_SORT_STATS_SAMTOOLS_SPLIT ( ch_split_corrected_bam,
                                        fasta )

        ch_split_sorted_bam = BAM_SORT_STATS_SAMTOOLS_SPLIT.out.bam
        ch_split_sorted_bai = BAM_SORT_STATS_SAMTOOLS_SPLIT.out.bai

        //
        // MODULE: Umitools Dedup
        //
        UMITOOLS_DEDUP ( ch_split_sorted_bam.join(ch_split_sorted_bai, by: [0]), true )

        ch_dedup_bam = UMITOOLS_DEDUP.out.bam
        ch_dedup_log = UMITOOLS_DEDUP.out.log
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        //
        // MODULE: Samtools merge
        //
        ch_bams_to_merge = ch_dedup_bam
                                .map{
                                    meta, bam ->
                                        bam_basename = bam.toString().split('/')[-1]
                                        split_bam_basename = bam_basename.split(/\./)
                                        meta = [ 'id': split_bam_basename[0] ]
                                    [ meta, bam ]
                                }
                                .groupTuple()

        SAMTOOLS_MERGE ( ch_bams_to_merge, fasta, fai)

        ch_dedup_merged_bam = SAMTOOLS_MERGE.out.bam

        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        BAM_SORT_STATS_SAMTOOLS_DEDUP ( ch_dedup_merged_bam,
                                        fasta )

        ch_dedup_sorted_bam = BAM_SORT_STATS_SAMTOOLS_DEDUP.out.bam
        ch_dedup_sorted_bai = BAM_SORT_STATS_SAMTOOLS_DEDUP.out.bai

        // these stats go for multiqc
        ch_dedup_sorted_stats = BAM_SORT_STATS_SAMTOOLS_DEDUP.out.stats
        ch_dedup_sorted_flagstat = BAM_SORT_STATS_SAMTOOLS_DEDUP.out.flagstat
        ch_dedup_sorted_idxstats = BAM_SORT_STATS_SAMTOOLS_DEDUP.out.idxstats
        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_DEDUP.out.versions)
    }

    //
    // MODULE: Isoquant
    //
    ISOQUANT ( ch_dedup_sorted_bam.join(ch_dedup_sorted_bai, by: [0]), gtf, fasta, fai, 'tag:CB')
    ch_gene_count_mtx = ISOQUANT.out.gene_count_mtx
    ch_transcript_count_mtx = ISOQUANT.out.transcript_count_mtx
    ch_versions = ch_versions.mix(ISOQUANT.out.versions)

    if (!params.skip_qc && !params.skip_seurat){
        //
        // MODULE: Seurat
        //
        SEURAT_GENE ( ch_gene_count_mtx.join(ch_dedup_sorted_flagstat, by: [0]) )
        ch_gene_seurat_qc = SEURAT_GENE.out.seurat_stats

        SEURAT_TRANSCRIPT ( ch_transcript_count_mtx.join(ch_dedup_sorted_flagstat, by: [0]) )
        ch_transcript_seurat_qc = SEURAT_TRANSCRIPT.out.seurat_stats

        //
        // MODULE: Combine Seurat Stats
        //

        ch_gene_stats = SEURAT_GENE.out.seurat_stats.collect{it[1]}
        COMBINE_SEURAT_STATS_GENE ( ch_gene_stats )
        ch_gene_stats_combined = COMBINE_SEURAT_STATS_GENE.out.combined_stats
        ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS_GENE.out.versions)

        ch_transcript_stats = SEURAT_TRANSCRIPT.out.seurat_stats.collect{it[1]}
        COMBINE_SEURAT_STATS_TRANSCRIPT ( ch_transcript_stats )
        ch_transcript_stats_combined = COMBINE_SEURAT_STATS_TRANSCRIPT.out.combined_stats
        ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS_TRANSCRIPT.out.versions)
    }

    //
    // SOFTWARE_VERSIONS
    //

    //
    // Collate and save software versions
    //
    //softwareVersionsToYAML(ch_versions)
    //    .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
    //    .set { ch_collated_versions }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_multiqc_report = Channel.empty()
    ch_versions = Channel.empty()

    if (!params.skip_qc && !params.skip_multiqc){

        //
        // MODULE: MultiQC for raw data
        //

        ch_multiqc_rawqc_files = Channel.empty()
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_config)
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
        // multiqc does not accept nanocomp with multiple samples as of now. See issue: https://github.com/ewels/MultiQC/issues/1630
        //ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_html.collect().ifEmpty([]))
        //ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.collect().ifEmpty([]))

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary    = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_finalqc_files = Channel.empty()
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_config)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_rseqc_read_dist.collect{it[1]}.ifEmpty([]))

        //ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_corrected_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_corrected_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_corrected_sorted_idxstats.collect{it[1]}.ifEmpty([]))

        //ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_sorted_idxstats.collect{it[1]}.ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_log.collect{it[1]}.ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_gene_stats_combined.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_transcript_stats_combined.collect().ifEmpty([]))

        MULTIQC_FINALQC (
            ch_multiqc_finalqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )
        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_versions    = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }

    emit:
    multiqc_report = ch_multiqc_report.toList()
    versions = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
