/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScnanoseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OTHER FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Dummy file used as optional input where required for __combine operator__
ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { NANOFILT                               } from "../modules/local/nanofilt"
include { PROWLERTRIMMER                         } from "../modules/local/prowlertrimmer"
include { SPLIT_FILE                             } from "../modules/local/split_file"
include { PIGZ as ZIP_R1                         } from "../modules/local/pigz"
include { PIGZ as ZIP_R2                         } from "../modules/local/pigz"
include { PIGZ as ZIP_TRIM                       } from "../modules/local/pigz"
include { PREEXTRACT_FASTQ                       } from "../modules/local/preextract_fastq"
include { UMI_TOOLS_WHITELIST                    } from "../modules/local/umi_tools_whitelist"
include { UMI_TOOLS_EXTRACT                      } from "../modules/local/umi_tools_extract"
include { PAFTOOLS                               } from "../modules/local/paftools"
include { MINIMAP2_INDEX                         } from "../modules/local/minimap2_index"
include { MINIMAP2_ALIGN                         } from "../modules/local/minimap2_align"
include { REFORMAT_WHITELIST                     } from "../modules/local/reformat_whitelist"
include { TAG_BARCODES                           } from "../modules/local/tag_barcodes"
include { CORRECT_BARCODES                       } from "../modules/local/correct_barcodes"
include { SORT_GTF                               } from "../modules/local/sort_gtf"
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPT } from '../modules/local/subread_featurecounts'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK             } from "../subworkflows/local/input_check"
include { CREATE_REGEX_INFO       } from "../subworkflows/local/create_regex"
include { PREPARE_REFERENCE_FILES } from "../subworkflows/local/prepare_reference_files"
include { GET_GENE_COUNTS_MATRIX  } from "../subworkflows/local/get_gene_counts_matrix"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP                                        } from "../modules/nf-core/gunzip/main"
include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { UMITOOLS_DEDUP                                } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BAM            } from "../modules/nf-core/samtools/view/main"
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER         } from "../modules/nf-core/samtools/view/main"
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP        } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BC_CORRECTED } from '../modules/nf-core/samtools/index/main'
include { STRINGTIE_STRINGTIE                           } from '../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE                               } from '../modules/nf-core/stringtie/merge/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM         } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM        } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_EXTRACTED    } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT     } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_MINIMAP  } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FILTERED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCNANOSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_fastq = INPUT_CHECK.out.reads

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - pre-trim QC
    //

    ch_fastqc_multiqc_pretrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_PRE_TRIM ( ch_fastq, params.skip_nanoplot, params.skip_fastqc)

        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // MODULE: Unzip fastq
    //
    GUNZIP( ch_fastq )
    ch_unzipped_fastqs = GUNZIP.out.gunzip

    //
    // SUBWORKFLOW: Prepare reference files
    //

    // TODO: Add fasta param to below
    PREPARE_REFERENCE_FILES ( "",
                            params.intron_retention_method,
                            params.fasta,
                            params.gtf)

    ch_fasta = PREPARE_REFERENCE_FILES.out.ch_prepared_fasta
    ch_gtf = PREPARE_REFERENCE_FILES.out.ch_prepared_gtf

    //
    // MODULE: Generate junction file - paftools
    //
    // TODO: *** once intron method 1/2 gets added, add conditionals to input gtf below (either param, or output of process) ***
    PAFTOOLS ( ch_gtf )
    ch_bed = PAFTOOLS.out.bed

    //
    // MODULE: Split fastq
    //
    ch_split_fastqs = ch_unzipped_fastqs

    if (params.split_amount > 0) {
        SPLIT_FILE( ch_unzipped_fastqs, '.fastq', params.split_amount)
        ch_split_fastqs = SPLIT_FILE.out.split_files

        // TODO: Change the ids so they contain the sample_name and index
    }

    ch_fastqs = Channel.empty()
    ch_split_fastqs.transpose().set { ch_fastqs }


    //
    // MODULE: Trim and filter reads
    //
    ch_trimmed_reads = ch_fastqs

    if (!params.skip_trimming){

        // TODO: Throw error if invalid trimming_software provided
        if (params.trimming_software == 'nanofilt') {

            NANOFILT ( ch_fastqs )
            ch_trimmed_reads = NANOFILT.out.reads

        } else if (params.trimming_software == 'prowler') {

            PROWLERTRIMMER ( ch_fastqs )
            ch_trimmed_reads = PROWLERTRIMMER.out.reads
        }
    }

    ch_trimmed_reads_qc = Channel.empty()
    ch_trimmed_reads
        .groupTuple()
        .set{ ch_trimmed_reads_qc }

    // TODO: may consider making zipped out below tmp. as it's only QC
    // NOTE: leaving it as is for now as we are testing
    ZIP_TRIM ( ch_trimmed_reads_qc, "filtered" )
    ch_zipped_trimmed_reads = ZIP_TRIM.out.archive

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
    //
    ch_fastqc_multiqc_postrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_POST_TRIM ( ch_zipped_trimmed_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // SUBWORKFLOW: Pre extract the cell barcodes
    //

    // We need to create the regex format
    // TODO: Add this information to the samplesheet to allow sample specific barcode detection?
    CREATE_REGEX_INFO( params.cell_barcode_pattern,
                params.identifier_pattern,
                params.cell_barcode_lengths,
                params.umi_lengths,
                params.fixed_seqs)

    val_regex_info = CREATE_REGEX_INFO.out.regex

    // Preextraction will create paired fastqs in cell ranger format
    // So we will need to set the fastqs to paired end
    ch_trimmed_reads
        .map {
            meta, fastq ->
                meta.single_end = false
                [ meta, fastq ]
        }

    PREEXTRACT_FASTQ( ch_trimmed_reads, val_regex_info.regex )

    // TODO: Cleaner way to do this and not repeat code?
    ch_pre_extracted_r1_fqs = Channel.empty()
    ch_pre_extracted_r2_fqs = Channel.empty()

    PREEXTRACT_FASTQ.out.r1_reads
        .groupTuple()
        .set{ ch_pre_extracted_r1_fqs }

    PREEXTRACT_FASTQ.out.r2_reads
        .groupTuple()
        .set{ ch_pre_extracted_r2_fqs }
    //
    // MODULE: Zip fastq
    //

    ZIP_R1 ( ch_pre_extracted_r1_fqs, "R1" )
    ch_zipped_r1_reads = ZIP_R1.out.archive

    ZIP_R2 ( ch_pre_extracted_r2_fqs, "R2" )
    ch_zipped_r2_reads = ZIP_R2.out.archive

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - pre-extracted QC
    //
    ch_fastqc_multiqc_pre_extracted = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_PRE_EXTRACTED ( ch_zipped_r2_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_pre_extracted = FASTQC_NANOPLOT_PRE_EXTRACTED.out.fastqc_multiqc.ifEmpty([])
    }

    // Merge the R1 and R2 fastqs back together
    ch_zipped_reads = Channel.empty()
    ch_zipped_r1_reads
        .join( ch_zipped_r2_reads )
        .map{ meta, r1, r2 ->
            [ meta, [r1, r2]]
        }
        .set{ ch_zipped_reads }

    //
    // MODULE: Create estimated whitelist
    //

    UMI_TOOLS_WHITELIST ( ch_zipped_reads, params.cell_amount, val_regex_info.umi_tools)
    ch_reads_with_whitelist = UMI_TOOLS_WHITELIST.out.whitelist

    //
    // MODULE: Extract barcodes
    //
    UMI_TOOLS_EXTRACT ( ch_reads_with_whitelist, val_regex_info.umi_tools )
    ch_extracted_reads = UMI_TOOLS_EXTRACT.out.reads

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_POST_EXTRACT ( ch_extracted_reads, params.skip_nanoplot, params.skip_fastqc)

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // MINIMAP2_INDEX
    //

    if (!params.skip_save_minimap2_index) {
        ch_fasta =  Channel.fromPath(params.fasta, checkIfExists: true)

        MINIMAP2_INDEX ( ch_fasta,  ch_bed)
        ch_minimap_index = MINIMAP2_INDEX.out.index
    }

    //
    // MINIMAP2_ALIGN
    //

    if (!params.skip_save_minimap2_index) {
        ch_reference = ch_minimap_index
    } else {
        ch_reference = Channel.fromPath(params.fasta, checkIfExists: true)
    }

    MINIMAP2_ALIGN ( ch_extracted_reads, ch_bed, ch_reference )

    MINIMAP2_ALIGN.out.sam
        .combine( ch_dummy_file )
        .set { ch_minimap_sam }

    //
    // MODULE: Samtools view
    //
    SAMTOOLS_VIEW_BAM ( ch_minimap_sam, [], [] )

    ch_minimap_bam = SAMTOOLS_VIEW_BAM.out.bam

    // acquire only mapped reads from bam for downstream processing

    ch_minimap_bam
        .combine( ch_dummy_file )
        .set { ch_minimap_bam_filter }

    SAMTOOLS_VIEW_FILTER ( ch_minimap_bam_filter, [], [] )
    ch_minimap_mapped_only_bam = SAMTOOLS_VIEW_FILTER.out.bam

    //
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
    BAM_SORT_STATS_SAMTOOLS_MINIMAP ( ch_minimap_bam, [] )
    ch_minimap_sorted_bam = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bam
    ch_minimap_sorted_bai = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bai
    // these stats go for multiqc
    ch_minimap_sorted_stats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.stats
    ch_minimap_sorted_flagstat = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.flagstat
    ch_minimap_sorted_idxstats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.idxstats

    BAM_SORT_STATS_SAMTOOLS_FILTERED ( ch_minimap_mapped_only_bam, [] )
    ch_minimap_filtered_sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam
    ch_minimap_filtered_sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai

    //
    // MODULE: Reformat whitelist
    //
    ch_whitelists = Channel.empty()
    ch_reads_with_whitelist
        .map{ meta, fastq, whitelist ->[ meta, whitelist ] }
        .set{ ch_whitelists }

    REFORMAT_WHITELIST ( ch_whitelists )
    ch_whitelist_bc_count = REFORMAT_WHITELIST.out.whitelist_bc_count
    ch_whitelist_bc_list = REFORMAT_WHITELIST.out.whitelist_bc_list

    //
    // MODULE: Tag Barcodes
    //

    ch_tag_barcode_in = Channel.empty()
    ch_minimap_filtered_sorted_bam
        .join( ch_zipped_r1_reads, by: 0 )
        .combine( val_regex_info.bc_length )
        .combine( val_regex_info.umi_length )
        .set{ ch_tag_barcode_in }

    TAG_BARCODES( ch_tag_barcode_in )
    ch_tagged_bam = TAG_BARCODES.out.tagged_bam

    //
    // MODULE: Correct Barcodes
    //

    ch_gt_whitelist = Channel.empty()
    if ( params.whitelist ) {
        ch_gt_whitelist = params.whitelist
    } else {
        ch_gt_whitelist = ch_whitelist_bc_list
    }

    ch_correct_barcode_in = Channel.empty()
    ch_tagged_bam
        .join ( ch_gt_whitelist, by: 0 )
        .join ( ch_whitelist_bc_count, by: 0 )
        .set { ch_correct_barcode_in }

    CORRECT_BARCODES( ch_correct_barcode_in )
    ch_corrected_bam = CORRECT_BARCODES.out.corrected_bam

    SAMTOOLS_INDEX_BC_CORRECTED ( ch_corrected_bam )
    ch_corrected_bam_bai = SAMTOOLS_INDEX_BC_CORRECTED.out.bai

    //
    // MODULE: Umitools Dedup
    //

    ch_dedup_bam = ch_corrected_bam
    ch_dedup_bam_bai = ch_corrected_bam_bai

    if (!params.skip_dedup) {
        UMITOOLS_DEDUP ( ch_corrected_bam.join(ch_corrected_bam_bai, by: [0]), true )

        ch_dedup_bam = UMITOOLS_DEDUP.out.bam

        SAMTOOLS_INDEX_DEDUP ( ch_dedup_bam )
        ch_dedup_bam_bai = SAMTOOLS_INDEX_DEDUP.out.bai
    }

    ch_dedup_bam
        .map{ meta, fastq ->
                meta.strandedness = params.stranded
                [ meta, fastq ]
            }
    //
    // SUBWORKFLOW: Gene Level Counts
    //

    // TODO: add 4th param (intron2 GTF) from PREPARE_REFERENCE_FILES
    if ( params.counts_level == 'gene' || !params.counts_level ) {
            if (!params.skip_dedup) {
                GET_GENE_COUNTS_MATRIX ( ch_dedup_bam, ch_gtf, params.intron_retention_method, "" )
                } else {
                GET_GENE_COUNTS_MATRIX ( ch_corrected_bam, ch_gtf, params.intron_retention_method, "" )
            }
        //TODO: change channed out to tagged matrix
        ch_gene_counts_mtx = GET_GENE_COUNTS_MATRIX.out.gene_counts_mtx
    }

    //TODO: bams below should be either dedub or not
    if ( params.counts_level == 'transcript' || !params.counts_level ) {
        //
        // MODULE: Create a stringtie gtf
        //
        STRINGTIE_STRINGTIE ( ch_dedup_bam, ch_gtf )
        ch_transcript_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf

        //
        // MODULE: Sort the gtf
        //
        SORT_GTF ( ch_transcript_gtf )
        ch_transcript_gtf_sorted = SORT_GTF.out.gtf

        //
        // MODULE: Merge the gtfs
        //
        // TODO: This currently doesn't take meta, so that means currently it does not work with multiple samples since files will get overwritten constantly. May want to convert this to local
        STRINGTIE_MERGE ( ch_transcript_gtf_sorted.map { meta, gtf -> gtf }, ch_gtf )
        ch_transcript_gtf_merged = STRINGTIE_MERGE.out.gtf

        //
        // MODULE: Transcript Level Counts
        //
        SUBREAD_FEATURECOUNTS_TRANSCRIPT ( ch_dedup_bam, ch_transcript_gtf_merged )
        ch_gene_transcript_mtx = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.counts
    }


    //
    // SOFTWARE_VERSIONS
    //

//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowScnanoseq.paramsSummaryMultiqc(workflow, summary_params)
//    ch_workflow_summary = Channel.value(workflow_summary)
//
//    ch_multiqc_files = Channel.empty()
//    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
//    //ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
//    //ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_pre_extracted.collect().ifEmpty([]))
//    //ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_minimap_sorted_stats.collect{it[1]}.ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_minimap_sorted_flagstat.collect{it[1]}.ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_minimap_sorted_idxstats.collect{it[1]}.ifEmpty([]))
//
//
//    MULTIQC (
//        ch_multiqc_files.collect()
//    )
//    multiqc_report = MULTIQC.out.report.toList()
//    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
