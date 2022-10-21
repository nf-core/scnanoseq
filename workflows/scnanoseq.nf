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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { NANOFILT            } from "../modules/local/nanofilt"
include { PROWLERTRIMMER      } from "../modules/local/prowlertrimmer"
include { SPLIT_FILE          } from "../modules/local/split_file"
include { PIGZ as ZIP_R1      } from "../modules/local/pigz"
include { PIGZ as ZIP_R2      } from "../modules/local/pigz"
include { PREEXTRACT_FASTQ    } from "../modules/local/preextract_fastq"
include { UMI_TOOLS_WHITELIST } from "../modules/local/umi_tools_whitelist"
include { UMI_TOOLS_EXTRACT } from "../modules/local/umi_tools_extract"
include { PAFTOOLS            } from "../modules/local/paftools.nf"
include { MINIMAP2_INDEX      } from "../modules/local/minimap2_index.nf"
include { MINIMAP2_ALIGN      } from "../modules/local/minimap2_align.nf"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK  } from '../subworkflows/local/input_check'
include { CREATE_REGEX_INFO } from "../subworkflows/local/create_regex"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP                      } from "../modules/nf-core/modules/gunzip/main"
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM     } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM    } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'

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

    // TODO: Turn trimming into subworkflow?
    //
    // MODULE: Unzip fastq
    //
    GUNZIP( ch_fastq )
    ch_unzipped_fastqs = GUNZIP.out.gunzip

    //
    // MODULE: Generate junction file - paftools
    //
    // TODO: *** once intron method 1/2 gets added, add conditionals to input gtf below (either param, or output of process) ***
    ch_gtf = file(params.gtf)
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

    //
    // SUBWORKFLOW: Pre extract the cell barcodes
    //

    // We need to create the regex format
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
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
    //
    ch_fastqc_multiqc_postrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_POST_TRIM ( ch_zipped_r2_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
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

    // TODO: How to read lines form the regex_pattern file?
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
    ch_minimap_index = MINIMAP2_ALIGN.out.sam
    ch_minimap_index.view()

    //
    // SOFTWARE_VERSIONS
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowScnanoseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
