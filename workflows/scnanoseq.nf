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

include { NANOFILT                                                 } from "../modules/local/nanofilt"
include { NANOCOMP as NANOCOMP_FASTQ                               } from "../modules/local/nanocomp"
include { NANOCOMP as NANOCOMP_BAM                                 } from "../modules/local/nanocomp"
include { PROWLERTRIMMER                                           } from "../modules/local/prowlertrimmer"
include { SPLIT_FILE                                               } from "../modules/local/split_file"
include { PIGZ as ZIP_R1                                           } from "../modules/local/pigz"
include { PIGZ as ZIP_R2                                           } from "../modules/local/pigz"
include { PIGZ as ZIP_TRIM                                         } from "../modules/local/pigz"
include { PREEXTRACT_FASTQ                                         } from "../modules/local/preextract_fastq"
include { UMI_TOOLS_WHITELIST                                      } from "../modules/local/umi_tools_whitelist"
include { UMI_TOOLS_EXTRACT                                        } from "../modules/local/umi_tools_extract"
include { PAFTOOLS                                                 } from "../modules/local/paftools"
include { MINIMAP2_INDEX                                           } from "../modules/local/minimap2_index"
include { MINIMAP2_ALIGN                                           } from "../modules/local/minimap2_align"
include { REFORMAT_WHITELIST                                       } from "../modules/local/reformat_whitelist"
include { TAG_BARCODES                                             } from "../modules/local/tag_barcodes"
include { CORRECT_BARCODES                                         } from "../modules/local/correct_barcodes"
include { SORT_GTF                                                 } from "../modules/local/sort_gtf"
include { MERGE_COUNTS_MTX                                         } from "../modules/local/merge_counts_mtx"
include { SEURAT as SEURAT_GENE                                    } from "../modules/local/seurat"
include { SEURAT as SEURAT_TRANSCRIPT                              } from "../modules/local/seurat"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_GENE        } from "../modules/local/combine_seurat_stats"
include { COMBINE_SEURAT_STATS as COMBINE_SEURAT_STATS_TRANSCRIPT  } from "../modules/local/combine_seurat_stats"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                           } from "../subworkflows/local/input_check"
include { CREATE_REGEX_INFO                                     } from "../subworkflows/local/create_regex"
include { PREPARE_REFERENCE_FILES                               } from "../subworkflows/local/prepare_reference_files"
include { GET_COUNTS_MATRIX as GET_GENE_COUNTS_MTX              } from "../subworkflows/local/get_counts_matrix"
include { GET_COUNTS_MATRIX as GET_TRANSCRIPT_COUNTS_MTX        } from "../subworkflows/local/get_counts_matrix"
include { GET_COUNTS_MATRIX as GET_INTRON_GENE_COUNTS_MTX       } from "../subworkflows/local/get_counts_matrix"
include { GET_COUNTS_MATRIX as GET_INTRON_TRANSCRIPT_COUNTS_MTX } from "../subworkflows/local/get_counts_matrix"

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

        FASTQC_NANOPLOT_PRE_TRIM ( ch_fastq, params.skip_nanoplot, params.skip_fastqc )

        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    if (!params.skip_qc && !params.skip_fastq_nanocomp) {
        ch_nanocomp_fastqs = ch_fastq.collect{it[1]}

        NANOCOMP_FASTQ ( ch_nanocomp_fastqs )
        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

    }

    //
    // MODULE: Unzip fastq
    //
    GUNZIP( ch_fastq )
    ch_unzipped_fastqs = GUNZIP.out.gunzip
    ch_versions = ch_versions.mix( GUNZIP.out.versions )

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES ( "",
                            params.intron_retention_method,
                            params.fasta,
                            params.gtf)

    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    ch_gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf
    ch_versions = ch_versions.mix( PREPARE_REFERENCE_FILES.out.versions )

    //
    // MODULE: Generate junction file - paftools
    //
    PAFTOOLS ( ch_gtf )
    ch_bed = PAFTOOLS.out.bed
    ch_versions = ch_versions.mix(PAFTOOLS.out.versions)

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
            .map{
                meta, fastq ->
                    def new_meta = [:]
                    new_meta.id = fastq.name.lastIndexOf('.').with {it != -1 ? fastq.name[0..<it] : fastq.name}
                    new_meta.single_end = meta.single_end
                [new_meta, fastq]
            }
            .set { ch_fastqs }

        ch_versions = ch_versions.mix(SPLIT_FILE.out.versions)
    }

    //
    // MODULE: Trim and filter reads
    //
    ch_trimmed_reads = ch_fastqs
    if (!params.skip_trimming){

        if (params.trimming_software == 'nanofilt') {

            NANOFILT ( ch_fastqs )
            ch_trimmed_reads = NANOFILT.out.reads
            ch_versions = ch_versions.mix(NANOFILT.out.versions)

        } else if (params.trimming_software == 'prowler') {

            PROWLERTRIMMER ( ch_fastqs )
            ch_trimmed_reads = PROWLERTRIMMER.out.reads
            ch_versions = ch_versions.mix(PROWLERTRIMMER.out.versions)

        }
        //
        // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
        //
        ch_fastqc_multiqc_postrim = Channel.empty()
        if (!params.skip_qc){
            // Concating the reads together temporarily for doing trim qc
            ch_trimmed_reads_qc = Channel.empty()
            ch_trimmed_reads
                .map {
                    meta, fastq ->
                        new_meta = [:]
                        new_meta.id = meta.id.split("\\.").with{it.length == 0 ? fastq.name : it[0..<-1].join('.')}
                        new_meta.single_end = true
                    [new_meta, fastq]
                }
                .groupTuple()
                .set { ch_trimmed_reads_qc }

            //
            // MODULE: Zip the reads
            //
            ZIP_TRIM (ch_trimmed_reads_qc, "filtered" )
            ch_zipped_trimmed_reads = ZIP_TRIM.out.archive
            ch_versions = ch_versions.mix(ZIP_TRIM.out.versions)

            //
            // MODULE: Run qc on the post trimmed reads
            //
            FASTQC_NANOPLOT_POST_TRIM ( ch_zipped_trimmed_reads, params.skip_nanoplot, params.skip_fastqc )

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.fastqc_version.first().ifEmpty(null))
        }
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
    // MODULE: Pre extract the cell barcodes
    //

    // Preextraction will create paired fastqs in cell ranger format
    // So we will need to set the fastqs to paired end

    PREEXTRACT_FASTQ( ch_trimmed_reads, val_regex_info.regex )

    ch_versions = ch_versions.mix(PREEXTRACT_FASTQ.out.versions)
    ch_pre_extracted_r1_fqs = Channel.empty() 
    ch_pre_extracted_r2_fqs = Channel.empty() 

    if ( params.split_amount > 0) {
        // TODO: Why does below work when the above solution doesn't?

        // THe reason we assign these variables is to ensure the order is consistent within the channels
        // TODO: Might be worth adding in a sorting to further make sure the sorting is correct
        ch_preextract_out_r1 = PREEXTRACT_FASTQ.out.r1_reads
        ch_preextract_out_r2 = PREEXTRACT_FASTQ.out.r2_reads

        ch_preextract_out_r1
            .map {
                meta, fastq ->
                    new_meta = ["id": meta.id.split("\\.").with{it.length == 0 ? fastq.name : it[0..<-1].join('.')},
                                "single_end": false]
                [new_meta, fastq]
            }
            .groupTuple()
            .set { ch_pre_extracted_r1_fqs }

        ch_preextract_out_r2
            .map {
                meta, fastq ->
                    new_meta = ["id": meta.id.split("\\.").with{it.length == 0 ? fastq.name : it[0..<-1].join('.')},
                                "single_end": false]
                [new_meta, fastq]
            }
            .groupTuple()
            .set { ch_pre_extracted_r2_fqs }

    } else {
        ch_pre_extracted_r1_fqs = PREEXTRACT_FASTQ.out.r1_reads 
        ch_pre_extracted_r2_fqs = PREEXTRACT_FASTQ.out.r2_reads

    }

    //
    // MODULE: Zip fastq
    //

    ZIP_R1 ( ch_pre_extracted_r1_fqs, "R1" )
    ch_zipped_r1_reads = ZIP_R1.out.archive
    ch_versions = ch_versions.mix(ZIP_R1.out.versions)

    ZIP_R2 ( ch_pre_extracted_r2_fqs, "R2" )
    ch_zipped_r2_reads = ZIP_R2.out.archive
    ch_versions = ch_versions.mix(ZIP_R2.out.versions)

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - pre-extracted QC
    //
    ch_fastqc_multiqc_pre_extracted = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_PRE_EXTRACTED ( ch_zipped_r2_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_pre_extracted = FASTQC_NANOPLOT_PRE_EXTRACTED.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_EXTRACTED.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_EXTRACTED.out.fastqc_version.first().ifEmpty(null))
    }

    // Merge the R1 and R2 fastqs back together
    ch_zipped_reads = Channel.empty()
    ch_zipped_r1_reads
        .join( ch_zipped_r2_reads )
        .map{ meta, r1, r2 ->
            [ meta, [r1, r2]]
        }
        .set{ ch_zipped_reads }

    
    ch_gt_whitelist = Channel.empty()
    if ( params.whitelist ) {
        ch_gt_whitelist = params.whitelist
    } else {
        //
        // MODULE: Create estimated whitelist
        //

        UMI_TOOLS_WHITELIST ( ch_zipped_reads, params.cell_amount, val_regex_info.umi_tools )
        ch_reads_with_whitelist = UMI_TOOLS_WHITELIST.out.whitelist
        ch_versions = ch_versions.mix(UMI_TOOLS_WHITELIST.out.versions)
        
        //
        // MODULE: Reformat whitelist
        //
        ch_whitelists = Channel.empty()
        ch_reads_with_whitelist
            .map{ meta, fastq, whitelist ->[ meta, whitelist ] }
            .set{ ch_whitelists }

        REFORMAT_WHITELIST ( ch_whitelists )
        ch_whitelist_bc_count = REFORMAT_WHITELIST.out.bc_list_counts
        ch_gt_whitelist = REFORMAT_WHITELIST.out.bc_list
        ch_versions = ch_versions.mix(REFORMAT_WHITELIST.out.versions)
        
    }

    //
    // MODULE: Extract barcodes
    //
    UMI_TOOLS_EXTRACT ( ch_zipped_reads.join(ch_gt_whitelist), val_regex_info.umi_tools )
    ch_extracted_reads = UMI_TOOLS_EXTRACT.out.reads
    ch_versions = ch_versions.mix(REFORMAT_WHITELIST.out.versions)

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_EXTRACT ( ch_extracted_reads, params.skip_nanoplot, params.skip_fastqc )

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_version.first().ifEmpty(null))
    }

    //
    // MINIMAP2_INDEX
    //

    if (!params.skip_save_minimap2_index) {
        ch_fasta =  Channel.fromPath(params.fasta, checkIfExists: true)

        MINIMAP2_INDEX ( ch_fasta,  ch_bed)
        ch_minimap_index = MINIMAP2_INDEX.out.index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MINIMAP2_ALIGN
    //

    if (!params.skip_save_minimap2_index) {
        ch_reference = ch_minimap_index.toList()
    } else {
        ch_reference = Channel.fromPath(params.fasta, checkIfExists: true).toList()
    }
    MINIMAP2_ALIGN ( ch_extracted_reads, ch_bed, ch_reference )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.sam
        .combine( ch_dummy_file )
        .set { ch_minimap_sam }

    //
    // MODULE: Samtools view
    //
    SAMTOOLS_VIEW_BAM ( ch_minimap_sam, [], [] )

    ch_minimap_bam = SAMTOOLS_VIEW_BAM.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_BAM.out.versions)

    // acquire only mapped reads from bam for downstream processing
    // NOTE: some QCs steps are performed on the full BAM

    ch_minimap_bam
        .combine( ch_dummy_file )
        .set { ch_minimap_bam_filter }

    SAMTOOLS_VIEW_FILTER ( ch_minimap_bam_filter, [], [] )
    ch_minimap_mapped_only_bam = SAMTOOLS_VIEW_FILTER.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILTER.out.versions)

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
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.versions)

    BAM_SORT_STATS_SAMTOOLS_FILTERED ( ch_minimap_mapped_only_bam, [] )
    ch_minimap_filtered_sorted_bam = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bam
    ch_minimap_filtered_sorted_bai = BAM_SORT_STATS_SAMTOOLS_FILTERED.out.bai
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_FILTERED.out.versions)

    //
    // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
    //

    if (!params.skip_qc && !params.skip_bam_nanocomp) {
        ch_nanocomp_bams = ch_minimap_sorted_bam.collect{it[1]}

        NANOCOMP_BAM ( ch_nanocomp_bams )
        ch_versions = ch_versions.mix( NANOCOMP_BAM.out.versions )
    }


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
    ch_versions = ch_versions.mix(TAG_BARCODES.out.versions)

    //
    // MODULE: Correct Barcodes
    //


    ch_correct_barcode_in = Channel.empty()
    ch_tagged_bam
        .join ( ch_gt_whitelist, by: 0 )
        .join ( ch_whitelist_bc_count, by: 0 )
        .set { ch_correct_barcode_in }

    CORRECT_BARCODES( ch_correct_barcode_in )
    ch_corrected_bam = CORRECT_BARCODES.out.corrected_bam
    ch_versions = ch_versions.mix(CORRECT_BARCODES.out.versions)

    SAMTOOLS_INDEX_BC_CORRECTED ( ch_corrected_bam )
    ch_corrected_bam_bai = SAMTOOLS_INDEX_BC_CORRECTED.out.bai
//    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_BC_CORRECTED)

    // TODO: Rename the dedup_bam channel to be more descriptive
    ch_dedup_bam = ch_corrected_bam
    ch_dedup_bam_bai = ch_corrected_bam_bai

    if (!params.skip_dedup) {
        //
        // MODULE: Umitools Dedup
        //
        UMITOOLS_DEDUP ( ch_corrected_bam.join(ch_corrected_bam_bai, by: [0]), true )
        ch_dedup_bam = UMITOOLS_DEDUP.out.bam
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        //
        // MODULE: Index the Dedup'd bam
        //
        SAMTOOLS_INDEX_DEDUP ( ch_dedup_bam )
        ch_dedup_bam_bai = SAMTOOLS_INDEX_DEDUP.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions)
    }

    ch_dedup_bam
        .map{ meta, fastq ->
                meta.strandedness = params.stranded
                [ meta, fastq ]
            }

    ch_gene_counts_mtx = Channel.empty()
    ch_transcript_counts_mtx = Channel.empty()

    if ( params.counts_level == "gene" || !params.counts_level ) {

        //
        // SUBWORKFLOW: Get the gene level count matrix
        //

        GET_GENE_COUNTS_MTX ( ch_dedup_bam, ch_gtf )
        ch_exon_gene_counts_mtx = GET_GENE_COUNTS_MTX.out.counts_mtx
        ch_gene_tag_bam_flagstat = GET_GENE_COUNTS_MTX.out.tag_bam_flagstat
        ch_versions = ch_versions.mix(GET_GENE_COUNTS_MTX.out.versions)

        if ( params.intron_retention_method == "2" ) {
            GET_INTRON_GENE_COUNTS_MTX ( ch_dedup_bam, ch_gtf )
            ch_intron_gene_counts_mtx = GET_INTRON_GENE_COUNTS_MTX.out.counts_mtx
            ch_versions = ch_versions.mix(GET_INTRON_GENE_COUNTS_MTX.out.versions)

            MERGE_COUNTS_MTX ( ch_exon_gene_counts_mtx.join ( ch_intron_gene_counts_mtx, by: 0 ))
            ch_gene_counts_mtx = MERGE_COUNTS_MTX.out.merged_mtx
            ch_versions = ch_versions.mix(MERGE_COUNTS_MTX.out.versions)
        } else {
            ch_gene_counts_mtx = ch_exon_gene_counts_mtx
        }
        //
        // MODULE: SEURAT
        //

        //TODO: may combine these into a single-channel later on

        ch_gene_stats_combined = Channel.empty()
        if (!params.skip_qc && !params.skip_seurat) {
            ch_gene_counts_flagstat = ch_gene_counts_mtx.join(ch_gene_tag_bam_flagstat, by: 0)
            SEURAT_GENE ( ch_gene_counts_flagstat )
            ch_versions = ch_versions.mix(SEURAT_GENE.out.versions)

            ch_gene_stats = SEURAT_GENE.out.seurat_stats.collect{it[1]}

            COMBINE_SEURAT_STATS_GENE ( ch_gene_stats )
            ch_gene_stats_combined = COMBINE_SEURAT_STATS_GENE.out.combined_stats
            ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS_GENE.out.versions)
        }

    }

    if ( params.counts_level == 'transcript' || !params.counts_level ) {
        //
        // MODULE: Create a stringtie gtf
        //
        STRINGTIE_STRINGTIE ( ch_dedup_bam, params.gtf )
        ch_transcript_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

        //
        // MODULE: Sort the gtf
        //
        SORT_GTF ( ch_transcript_gtf )
        ch_transcript_gtf_sorted = SORT_GTF.out.gtf
        ch_versions = ch_versions.mix(SORT_GTF.out.versions)

        //
        // MODULE: Merge the gtfs
        //
        // TODO: This currently doesn't take meta, so that means currently it does not work with multiple samples since files will get overwritten constantly. May want to convert this to local
        STRINGTIE_MERGE ( ch_transcript_gtf_sorted.map { meta, gtf -> gtf }, ch_gtf )
        ch_transcript_gtf_merged = STRINGTIE_MERGE.out.gtf
        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)

        //
        // SUBWORKFLOW: Get the transcript level count matrix
        //
        GET_TRANSCRIPT_COUNTS_MTX ( ch_dedup_bam, ch_transcript_gtf_merged )
        ch_transcript_counts_mtx = GET_TRANSCRIPT_COUNTS_MTX.out.counts_mtx
        ch_transcript_tag_bam_flagstat = GET_TRANSCRIPT_COUNTS_MTX.out.tag_bam_flagstat
        ch_versions = ch_versions.mix(GET_TRANSCRIPT_COUNTS_MTX.out.versions)
        //
        // MODULE: SEURAT
        //

        //TODO: may combine these into a single-channel later on

        ch_transcript_stats_combined = Channel.empty()
        if (!params.skip_qc && !params.skip_seurat) {
            ch_transcript_counts_flagstat = ch_transcript_counts_mtx.join(ch_transcript_tag_bam_flagstat, by: 0)

            SEURAT_TRANSCRIPT ( ch_transcript_counts_flagstat )

            ch_transcript_stats = SEURAT_TRANSCRIPT.out.seurat_stats.collect{it[1]}
            ch_versions = ch_versions.mix(SEURAT_TRANSCRIPT.out.versions)

            COMBINE_SEURAT_STATS_TRANSCRIPT ( ch_transcript_stats )
            ch_transcript_stats_combined = COMBINE_SEURAT_STATS_TRANSCRIPT.out.combined_stats
            ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS_TRANSCRIPT.out.versions)
        }

    }


    // TODO: combine seurat stats to generat MultiQC table

    //
    // SOFTWARE_VERSIONS
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (!params.skip_multiqc){

        //
        // MODULE: MultiQC for raw data
        //

        ch_multiqc_rawqc_files = Channel.empty()
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect()
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        workflow_summary    = WorkflowScnanoseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_finalqc_files = Channel.empty()
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_gene_stats_combined.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_transcript_stats_combined.collect().ifEmpty([]))

        MULTIQC_FINALQC (
            ch_multiqc_finalqc_files.collect()
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
