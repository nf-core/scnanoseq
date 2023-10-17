/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowScnanoseq.initialise(params, log)

def checkPathParamList = [ 
    params.input, params.multiqc_config, params.fasta, 
    params.gtf, params.whitelist
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER PRESETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This is if the user passes in direct regex
def cell_barcode_pattern = ""

// This is for if the user wants to do more human readable regex
def identifier_pattern = ""
def cell_barcode_lengths = ""
def umi_lengths = ""
def fixed_seqs = ""

if (params.barcode_preset) {
    if (params.barcode_preset = "cellranger_3_prime") {
        identifier_pattern = "fixed_seq_1,cell_barcode_1,umi_1,fixed_seq_2"
        cell_barcode_lengths = "16"
        umi_lengths = "12"
        fixed_seqs = "CTACACGACGCTCTTCCGATCT, TTTTTTTTTT"

    } else if (params.barcode_preset = "cellranger_5_prime") {
        identifier_pattern = "fixed_seq_1,cell_barcode_1,umi_1,fixed_seq_2"
        cell_barcode_lengths = "16"
        umi_lengths = "12"
        fixed_seqs = "CTACACGACGCTCTTCCGATCT, TTTTTTTTTT"
    }
} else {
    identifier_pattern = params.identifier_pattern
    cell_barcode_lengths = params.cell_barcode_lengths
    umi_lengths = params.umi_lengths
    fixed_seqs = params.fixed_seqs

}

// TODO: Adding this in temporarily. Rethink how we want to represent this
def blaze_whitelist = file("$baseDir/assets/whitelist/3M-february-2018.zip")
if (params.whitelist) {
    blaze_whitelist = whitelist
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
include { PROWLERTRIMMER                                           } from "../modules/local/prowlertrimmer"
include { SPLIT_FILE                                               } from "../modules/local/split_file"
include { PIGZ as ZIP_TRIM                                         } from "../modules/local/pigz"
include { BLAZE                                                    } from "../modules/local/blaze"
include { PREEXTRACT_FASTQ                                         } from "../modules/local/preextract_fastq.nf"
include { PAFTOOLS                                                 } from "../modules/local/paftools"
include { MINIMAP2_INDEX                                           } from "../modules/local/minimap2_index"
include { MINIMAP2_ALIGN                                           } from "../modules/local/minimap2_align"
include { TAG_BARCODES                                             } from "../modules/local/tag_barcodes"
include { CORRECT_BARCODES                                         } from "../modules/local/correct_barcodes"
include { MERGE_COUNTS_MTX                                         } from "../modules/local/merge_counts_mtx"
include { ISOQUANT                                                 } from "../modules/local/isoquant"
include { SEURAT as SEURAT_GENE                                    } from "../modules/local/seurat"
include { SEURAT as SEURAT_TRANSCRIPT                              } from "../modules/local/seurat"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_GENE        } from "../modules/local/combine_seurat_stats"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_TRANSCRIPT  } from "../modules/local/combine_seurat_stats"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK             } from "../subworkflows/local/input_check"
include { CREATE_REGEX_INFO       } from "../subworkflows/local/create_regex"
include { PREPARE_REFERENCE_FILES } from "../subworkflows/local/prepare_reference_files"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP                                        } from "../modules/nf-core/gunzip/main"
include { MULTIQC as MULTIQC_RAWQC                      } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_FINALQC                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { UMITOOLS_DEDUP                                } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BAM            } from "../modules/nf-core/samtools/view/main"
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER         } from "../modules/nf-core/samtools/view/main"
include { CAT_CAT                                       } from "../modules/nf-core/cat/cat/main"
include { CAT_FASTQ                                     } from '../modules/nf-core/cat/fastq/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM         } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM        } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT     } from '../subworkflows/nf-core/qcfastq_nanoplot_fastqc'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_MINIMAP  } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FILTERED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_CORRECTED } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_DEDUP } from "../subworkflows/nf-core/bam_sort_stats_samtools/main"

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
    INPUT_CHECK ( ch_input )
        .reads
        .map{
            meta, fastq ->
                new_id = meta.id - ~/_T\d+/
                [ meta + [id: new_id], fastq ]
        }
        .groupTuple()
        .branch{
            meta, fastq ->
                single: fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastqs }

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    //
    // MODULE: Combine fastqs from the same sample
    //
    CAT_FASTQ ( ch_fastqs.multiple )
        .reads
        .mix ( ch_fastqs.single )
        .set { ch_cat_fastq }

    ch_versions = ch_versions.mix (CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - pre-trim QC
    //

    ch_fastqc_multiqc_pretrim = Channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_PRE_TRIM ( ch_cat_fastq, params.skip_nanoplot, params.skip_fastqc )

        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = Channel.empty()
    ch_nanocomp_fastq_txt = Channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {
        ch_nanocomp_fastqs = ch_cat_fastq.collect{it[1]}

        NANOCOMP_FASTQ ( ch_nanocomp_fastqs )
        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.txt

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

    }

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES ( "",
                            params.intron_retention_method,
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

        if (params.trimming_software == 'nanofilt') {

            NANOFILT ( ch_fastqs )
            ch_trimmed_reads = NANOFILT.out.reads
            ch_versions = ch_versions.mix(NANOFILT.out.versions)

        } else if (params.trimming_software == 'prowler') {

            PROWLERTRIMMER ( ch_fastqs )
            ch_trimmed_reads = PROWLERTRIMMER.out.reads
            ch_versions = ch_versions.mix(PROWLERTRIMMER.out.versions)

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
            FASTQC_NANOPLOT_POST_TRIM ( ch_zipped_reads, params.skip_nanoplot, params.skip_fastqc )

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.fastqc_version.first().ifEmpty(null))
        }
    } else {
        ch_zipped_reads = ch_cat_fastq 
    }
    
    //
    // MODULE: Parse the regex info
    //

    // We need to create the regex format
    // TODO: Add this information to the samplesheet to allow sample specific barcode detection?
    CREATE_REGEX_INFO( cell_barcode_pattern,
                identifier_pattern,
                cell_barcode_lengths,
                umi_lengths,
                fixed_seqs)

    val_regex_info = CREATE_REGEX_INFO.out.regex
    // TODO: Why can't we use the below code?
    //ch_versions = ch_versions.mix(CREATE_REGEX_INFO.out.versions)

    //
    // MODULE: Generate whitelist
    //

    BLAZE ( ch_zipped_reads, params.cell_amount, blaze_whitelist)

    ch_putative_bc = BLAZE.out.putative_bc
    ch_gt_whitelist = BLAZE.out.whitelist
    ch_whitelist_bc_count = BLAZE.out.bc_count
    ch_versions = ch_versions.mix(BLAZE.out.versions)

    //
    // MODULE: Extract barcodes
    //

    PREEXTRACT_FASTQ( ch_zipped_reads.join(ch_putative_bc))
    ch_zipped_r1_reads = PREEXTRACT_FASTQ.out.r1_reads
    ch_zipped_r2_reads = PREEXTRACT_FASTQ.out.r2_reads

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_EXTRACT ( ch_zipped_r2_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_version.first().ifEmpty(null))
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
    // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
    //
    ch_nanocomp_bam_html = Channel.empty()
    ch_nanocomp_bam_txt = Channel.empty()

    if (!params.skip_qc && !params.skip_bam_nanocomp) {
        ch_nanocomp_bams = ch_minimap_sorted_bam.collect{it[1]}

        NANOCOMP_BAM ( ch_nanocomp_bams )
        
        ch_nanocomp_bam_html = NANOCOMP_BAM.out.html 
        ch_nanocomp_bam_txt = NANOCOMP_BAM.out.txt
        ch_versions = ch_versions.mix( NANOCOMP_BAM.out.versions )
    }


    //
    // MODULE: Tag Barcodes
    //

    TAG_BARCODES (
        ch_minimap_filtered_sorted_bam
            .join( ch_zipped_r1_reads, by: 0 )
            .combine( val_regex_info.bc_length )
            .combine( val_regex_info.umi_length )
    )

    ch_tagged_bam = TAG_BARCODES.out.tagged_bam
    ch_versions = ch_versions.mix(TAG_BARCODES.out.versions)

    //
    // MODULE: Correct Barcodes
    //

    CORRECT_BARCODES (
        ch_tagged_bam
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
        // MODULE: Umitools Dedup
        //
        UMITOOLS_DEDUP ( ch_corrected_sorted_bam.join(ch_corrected_sorted_bai, by: [0]), true )

        ch_dedup_bam = UMITOOLS_DEDUP.out.bam
        ch_dedup_log = UMITOOLS_DEDUP.out.log
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
        // The subworkflow is called in both the minimap2 bams and filtered (mapped only) version
        BAM_SORT_STATS_SAMTOOLS_DEDUP ( ch_dedup_bam,
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

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (!params.skip_qc && !params.skip_multiqc){

        //
        // MODULE: MultiQC for raw data
        //

        ch_multiqc_rawqc_files = Channel.empty()
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_config)
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_html.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.collect().ifEmpty([]))

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        workflow_summary    = WorkflowScnanoseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_finalqc_files = Channel.empty()
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_config)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_nanocomp_bam_html.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_nanocomp_bam_txt.collect().ifEmpty([]))
        
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_idxstats.collect{it[1]}.ifEmpty([]))
        
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
        multiqc_report = MULTIQC_FINALQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
