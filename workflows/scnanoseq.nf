/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER PRESETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Whitelist
if (params.cdna_whitelist) {
    cdna_whitelist = file(params.cdna_whitelist)
}
else {
    if (params.barcode_format.equals("10X_3v3")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-february-2018.txt.gz")
    }
    else if (params.barcode_format.equals("10X_5v2")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/737K-august-2016.txt.gz")
    }
    else if (params.barcode_format.equals("10X_3v4")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-3pgex-may-2023_TRU.txt.gz")
    }
    else if (params.barcode_format.equals("10X_5v3")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-5pgex-jan-2023.txt.gz")
    }
    else if (params.barcode_format.equals("10X_multiome")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/cellranger_arc_rna.737K-arc-v1.txt.gz")
    } else {
        cdna_whitelist = []
    }
}

if (params.dna_whitelist) {
    dna_whitelist = file(params.dna_whitelist)
}
else if (params.barcode_format.equals("10X_multiome")) {
    dna_whitelist = file("$baseDir/assets/whitelist/cellranger_arc_atac.737K-arc-v1.txt.gz")
} else {
    dna_whitelist = []
}

// Quantifiers

// Associate the quantifiers with the kind of alignment needed
GENOME_QUANT_OPTS = [ 'isoquant' ]
TRANSCRIPT_QUANT_OPTS = [ 'oarfish' ]

genome_quants = []
transcript_quants = []
for (quantifier in params.quantifier.split(',')) {
    if (quantifier in GENOME_QUANT_OPTS) {
        genome_quants.add(quantifier)
    }

    if (quantifier in TRANSCRIPT_QUANT_OPTS) {
        transcript_quants.add(quantifier)
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                       = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config                = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                         = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description   = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { NANOFILT                          } from "../modules/local/nanofilt"
include { SPLIT_FILE                        } from "../modules/local/split_file"
include { PREEXTRACT_FASTQ                  } from "../modules/local/preextract_fastq.nf"
include { READ_COUNTS                       } from "../modules/local/read_counts.nf"
include { CORRECT_BARCODES                  } from "../modules/local/correct_barcodes"
include { UCSC_GTFTOGENEPRED                } from "../modules/local/ucsc_gtftogenepred"
include { UCSC_GENEPREDTOBED                } from "../modules/local/ucsc_genepredtobed"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PREPARE_REFERENCE_FILES                                     } from "../subworkflows/local/prepare_reference_files"
include { DEMULTIPLEX_FLEXIPLEX as DEMULTIPLEX_FLEXIPLEX_CDNA         } from "../subworkflows/local/demultiplex_flexiplex"
include { DEMULTIPLEX_FLEXIPLEX as DEMULTIPLEX_FLEXIPLEX_DNA          } from "../subworkflows/local/demultiplex_flexiplex"
include { DEMULTIPLEX_BLAZE                                           } from "../subworkflows/local/demultiplex_blaze"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_GENOME     } from "../subworkflows/local/process_longread_scrna"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_TRANSCRIPT } from "../subworkflows/local/process_longread_scrna"
include { ALIGN_DEDUPLICATE_DNA                                       } from "../subworkflows/local/align_deduplicate_dna"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//

include { NANOCOMP as NANOCOMP_FASTQ                    } from "../modules/nf-core/nanocomp/main"
include { CHOPPER                                       } from "../modules/nf-core/chopper/main"
include { MULTIQC as MULTIQC_RAWQC                      } from "../modules/nf-core/multiqc/main"
include { MULTIQC as MULTIQC_FINALQC                    } from "../modules/nf-core/multiqc/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from "../modules/nf-core/custom/dumpsoftwareversions/main"
include { CAT_CAT                                       } from "../modules/nf-core/cat/cat/main"
include { CAT_FASTQ                                     } from "../modules/nf-core/cat/fastq/main"
include { paramsSummaryMap                              } from "plugin/nf-schema"

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/subworkflows
 */

include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM          } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM         } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT      } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { paramsSummaryMultiqc                                         } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML                                       } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText                                       } from "../subworkflows/local/utils_nfcore_scnanoseq_pipeline"

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
    ch_multiqc_report = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    ch_samplesheet
        .branch {
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
        ch_nanostat_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_txt.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = Channel.empty()
    ch_nanocomp_fastq_txt = Channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ (
            ch_cat_fastq
                .collect{it[1]}
                .map{
                    [ [ 'id': 'nanocomp_fastq.' ] , it ]
                }
        )

        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.report_html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.stats_txt

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

    }

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES (
        params.genome_fasta,
        params.transcript_fasta,
        params.gtf
    )

    genome_fasta = PREPARE_REFERENCE_FILES.out.prepped_genome_fasta
    genome_fai = PREPARE_REFERENCE_FILES.out.genome_fai
    transcript_fasta = PREPARE_REFERENCE_FILES.out.prepped_transcript_fasta
    transcript_fai = PREPARE_REFERENCE_FILES.out.transcript_fai
    gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf

    ch_versions = ch_versions.mix( PREPARE_REFERENCE_FILES.out.versions )

    //
    // MODULE: Generate bed file from input gtf for rseqc
    //

    // come back to this once intron work is finished (likely input will be fine)
    ch_pred = Channel.empty()
    ch_rseqc_bed = Channel.empty()
    if (!params.skip_qc && !params.skip_rseqc) {
        UCSC_GTFTOGENEPRED( gtf )
        ch_pred = UCSC_GTFTOGENEPRED.out.genepred
        ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

        UCSC_GENEPREDTOBED ( ch_pred )
        ch_rseqc_bed = UCSC_GENEPREDTOBED.out.bed
        ch_versions = ch_versions.mix(UCSC_GENEPREDTOBED.out.versions)
    }

    //
    // MODULE: Trim and filter reads
    //
    ch_fastqc_multiqc_postrim = Channel.empty()
    ch_trimmed_reads_combined = Channel.empty()

    if (!params.skip_trimming){

        //
        // MODULE: Chopper
        //

        CHOPPER ( ch_cat_fastq, [] )

        versions = CHOPPER.out.versions
        ch_trimmed_reads_combined = CHOPPER.out.fastq

        //
        // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
        //
        if (!params.skip_qc){

            //
            // MODULE: Run qc on the post trimmed reads
            //
            FASTQC_NANOPLOT_POST_TRIM ( ch_trimmed_reads_combined, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_nanostat_posttrim = FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_txt.ifEmpty([])
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.toulligqc_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.fastqc_version.first().ifEmpty(null))
        }
    } else {
        ch_trimmed_reads_combined = ch_unzipped_fastqs
    }


    // Branch channel to dna and cdna
    ch_trimmed_reads_combined = ch_trimmed_reads_combined
        .branch {
            meta, fastq ->
                dna: meta.type == 'dna'
                    return [ meta, fastq ]
                cdna: meta.type == 'cdna'
                    return [ meta, fastq ]
        }

    //
    // SUBWORKFLOW: Demultiplex reads using FLEXIPLEX for DNA
    //

    //TODO: Add check for nonempty dna channels not whitelist presence

    DEMULTIPLEX_FLEXIPLEX_DNA (
        ch_trimmed_reads_combined.dna,
        dna_whitelist
    )

    ch_versions = ch_versions.mix(DEMULTIPLEX_FLEXIPLEX_DNA.out.versions)
    ch_extracted_fastq_dna = DEMULTIPLEX_FLEXIPLEX_DNA.out.flexiplex_fastq
    ch_corrected_bc_info_dna = DEMULTIPLEX_FLEXIPLEX_DNA.out.flexiplex_barcodes


    if (params.demux_tool == "flexiplex") {

        //
        // SUBWORKFLOW: Demultiplex reads using FLEXIPLEX for cDNA
        //

        DEMULTIPLEX_FLEXIPLEX_CDNA (
            ch_trimmed_reads_combined.cdna,
            cdna_whitelist
        )

        ch_versions = ch_versions.mix(DEMULTIPLEX_FLEXIPLEX_CDNA.out.versions)

        ch_extracted_fastq_cdna = DEMULTIPLEX_FLEXIPLEX_CDNA.out.flexiplex_fastq
        ch_corrected_bc_info_cdna = DEMULTIPLEX_FLEXIPLEX_CDNA.out.flexiplex_barcodes

    } else if (params.demux_tool == "blaze") {

        //
        // SUBWORKFLOW: Demultiplex reads using BLAZE for cDNA
        //

        DEMULTIPLEX_BLAZE (
            ch_trimmed_reads_combined.cdna,
            cdna_whitelist
        )

        ch_versions = ch_versions.mix(DEMULTIPLEX_BLAZE.out.versions)

        ch_extracted_fastq_cdna = DEMULTIPLEX_BLAZE.out.extracted_fastq
        ch_corrected_bc_info_cdna = DEMULTIPLEX_BLAZE.out.corrected_bc_info
    }

    // Recombine channels for QC modules
    ch_extracted_fastq = ch_extracted_fastq_cdna.mix(ch_extracted_fastq_dna)
    ch_corrected_bc_info = ch_corrected_bc_info_cdna.mix(ch_corrected_bc_info_dna)

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = Channel.empty()
    ch_read_counts = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_EXTRACT ( ch_extracted_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_nanostat_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_txt.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_version.first().ifEmpty(null))

        //
        // MODULE: Generate read counts
        //

        ch_pretrim_counts = Channel.empty()
        ch_posttrim_counts = Channel.empty()
        ch_postextract_counts = Channel.empty()
        if (!params.skip_fastqc){
            ch_pretrim_counts = ch_fastqc_multiqc_pretrim.collect{it[0]}
            ch_posttrim_counts = ch_fastqc_multiqc_postrim.collect{it[0]}
            ch_postextract_counts = ch_fastqc_multiqc_postextract.collect{it[0]}

        } else if (!params.skip_nanoplot){
            ch_pretrim_counts = ch_nanostat_pretrim.collect{it[1]}
            ch_posttrim_counts = ch_nanostat_posttrim.collect{it[1]}
            ch_postextract_counts = ch_nanostat_postextract.collect{it[1]}

        }

        READ_COUNTS (
            ch_pretrim_counts.ifEmpty([]),
            ch_posttrim_counts.ifEmpty([]),
            ch_postextract_counts.ifEmpty([]),
            ch_corrected_bc_info.collect{it[1]})

        ch_read_counts = READ_COUNTS.out.read_counts
        ch_versions = ch_versions.mix(READ_COUNTS.out.versions)
    }

    //
    // SUBWORKFLOW: Align Long Read cDNA data
    //

    ch_multiqc_finalqc_files = Channel.empty()

    if (genome_quants){
        PROCESS_LONGREAD_SCRNA_GENOME(
            genome_fasta,
            genome_fai,
            gtf,
            ch_extracted_fastq_cdna,
            ch_rseqc_bed,
            ch_corrected_bc_info_cdna,
            genome_quants,
            params.dedup_tool,
            true, // Used to indicate the bam is genome aligned
            params.fasta_delimiter,
            params.skip_save_minimap2_index,
            params.skip_qc,
            params.skip_rseqc,
            params.skip_bam_nanocomp,
            params.skip_seurat,
            params.skip_dedup
        )
        ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_GENOME.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_flagstat.collect{it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_idxstats.collect{it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_rseqc_read_dist.collect{it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_nanocomp_bam_txt.collect{it[1]}.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.bc_tagged_flagstat.collect{it[1]}.ifEmpty([])
        )

        if (!params.skip_dedup) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_flagstat.collect{it[1]}.ifEmpty([])
            )
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_idxstats.collect{it[1]}.ifEmpty([])
            )
        }

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            ch_read_counts.collect().ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.gene_qc_stats.collect().ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.transcript_qc_stats.collect().ifEmpty([])
        )
    }

    // oarfish expects deduplicated reads
    if (transcript_quants) {
        PROCESS_LONGREAD_SCRNA_TRANSCRIPT (
            transcript_fasta,
            transcript_fai,
            gtf,
            ch_extracted_fastq_cdna,
            ch_rseqc_bed,
            ch_corrected_bc_info_cdna,
            transcript_quants,
            params.dedup_tool,
            false, // Indicates this is NOT genome aligned
            params.fasta_delimiter,
            params.skip_save_minimap2_index,
            params.skip_qc,
            true, // RSeQC does not work well with transcriptome alignments
            true, // Nanocomp does not work well with transcriptome alignments
            params.skip_seurat,
            false // Oarfish requires deduplication, so cannot skip it
        )

        ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_flagstat.collect{it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_rseqc_read_dist.collect{it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_nanocomp_bam_txt.collect{it[1]}.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.bc_tagged_flagstat.collect{it[1]}.ifEmpty([])
        )

        if (!params.skip_dedup) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.dedup_flagstat.collect{it[1]}.ifEmpty([])
            )
        }

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            ch_read_counts.collect().ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.transcript_qc_stats.collect().ifEmpty([])
        )
    }

    //
    // SUBWORKFLOW: Align and deduplicate DNA samples
    //

    ALIGN_DEDUPLICATE_DNA (
        genome_fasta,
        genome_fai,
        ch_extracted_fastq_dna,
        params.skip_save_minimap2_index,
        params.skip_qc,
        params.skip_bam_nanocomp,
        params.skip_dedup
    )

    ch_versions = ch_versions.mix(ALIGN_DEDUPLICATE_DNA.out.versions)
    
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.flagstat.collect{it[1]}.ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.nanocomp_bam_txt.collect{it[1]}.ifEmpty([])
    )

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
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.collect{it[1]}.ifEmpty([]))

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([]),
            [],
            []
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary    = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_config)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))


        MULTIQC_FINALQC (
            ch_multiqc_finalqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([]),
            [],
            []
        )
        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_versions    = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }

    emit:
    multiqc_report = ch_multiqc_report.toList()
    versions       = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
