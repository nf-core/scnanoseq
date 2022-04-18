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

include { NANOFILT                   } from "../modules/local/nanofilt"
include { PROWLERTRIMMER             } from "../modules/local/prowlertrimmer"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM        } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM       } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT    } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'

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
    // MODULE: Trim and filter reads
    //

    // TODO: Throw error if invalid trimming_software provided
    ch_trimmed_reads = ch_fastq
    if (!params.skip_trimming){

        if (params.trimming_software == 'nanofilt') {

            NANOFILT ( ch_fastq )
            ch_trimmed_reads = NANOFILT.out.reads
        } else if (params.trimming_software == 'prowler') {

            PROWLERTRIMMER ( ch_fastq )
            ch_trimmed_reads = PROWLERTRIMMER.out.reads
        }
    }

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
    //
    //TODO: enable below for post-trim QC; change input channel <<<ch_fastq>>>
    /*
    ch_fastqc_multiqc_postrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_POST_TRIM ( <<<ch_fastq>>>, params.skip_nanoplot, params.skip_fastqc)

        ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
    }
    */

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    //TODO: enable below for post-extract QC;  change input channel <<<ch_fastq>>>
    /*
    ch_fastqc_multiqc_postextract = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_POST_EXTRACT ( <<<ch_fastq>>>, params.skip_nanoplot, params.skip_fastqc)

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
    }
    */

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
