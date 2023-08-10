//
// Creates gtfs to that add introns as features
//


// Local modules
include { TRANSCRIPT_TO_EXON        } from '../../modules/local/prepare_gtf_transcript_to_exon'
include { GET_GTF_FEATURES          } from '../../modules/local/get_gtf_features'
include { GTF2BED                   } from '../../modules/local/gtf2bed'
include { UCSC_BEDTOGENEPRED        } from '../../modules/local/ucsc_bedtogenepred'
include { UCSC_GENEPREDTOGTF        } from '../../modules/local/ucsc_genepredtogtf'
include { CREATE_INTRON_GTF         } from '../../modules/local/create_intron_gtf'

// nf-core modules
include { CUSTOM_GETCHROMSIZES                        } from '../../modules/nf-core/custom/getchromsizes/main'
include { BEDTOOLS_COMPLEMENT as COMPLEMENT_GTF       } from '../../modules/nf-core/bedtools/complement/main'
include { BEDTOOLS_COMPLEMENT as COMPLEMENT_NONINTRON } from '../../modules/nf-core/bedtools/complement/main'
include { CAT_CAT as CAT_BED                          } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_GTF                          } from '../../modules/nf-core/cat/cat/main'

workflow PREPARE_GTF {
    take:
    gtf_preparation_method
    gtf
    fasta

    main:

    ch_versions = Channel.empty()
    ch_prepared_gtf = Channel.empty()

    if (gtf_preparation_method == "1") {
        // Convert the Transcript to exons
        TRANSCRIPT_TO_EXON(gtf)
        ch_prepared_gtf = TRANSCRIPT_TO_EXON.out.ch_processed_gtf
        ch_versions = ch_versions.mix(TRANSCRIPT_TO_EXON.out.versions)

    } else if (gtf_preparation_method == "2") {
        //
        // MODULE: Get the chromosome sizes
        //
        CUSTOM_GETCHROMSIZES ( [ ["id":"chr_sizes"], fasta ])
        chr_sizes = CUSTOM_GETCHROMSIZES.out.sizes
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

        //
        // MODULE: Complement the gtf to get the intergenic regions
        //

        COMPLEMENT_GTF (
            gtf
                .map {
                    meta, gtf ->
                        meta.id = "intergenic"
                        [meta, gtf]
                },
            chr_sizes 
        )
        ch_intergenic_bed = COMPLEMENT_GTF.out.bed
        ch_versions = ch_versions.mix(COMPLEMENT_GTF.out.versions)

        //
        // MODULE: Get the exon regions
        //
        GET_GTF_FEATURES( [ [ "id": "exon" ], gtf], "exon" )
        ch_exon_gtf = GET_GTF_FEATURES.out.gtf
        ch_versions = ch_versions.mix(GET_GTF_FEATURES.out.versions)

        //
        // MODULE: Convert the gtf to bed
        //
        GTF2BED( ch_exon_gtf )
        ch_exon_bed = GTF2BED.out.bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)

        // Need to format the file list so that they are a single list
        ch_files = Channel.empty()
        ch_intergenic_bed
            .map { meta, bed -> [bed] }
            .concat (
                ch_exon_bed
                    .map { meta, bed -> [bed] }
                )
            .collect()
            .toList()
            .set{ch_files}

        //
        // MODULE: Combine the file list with a meta object
        //
        ch_cat_files_in = Channel.of(["id": "not_introns"]).concat(ch_files).collect()

        CAT_BED( ch_cat_files_in )

        ch_not_intron_bed = CAT_BED.out.file_out
        ch_versions = ch_versions.mix(CAT_BED.out.versions)

        //
        // MODULE: Complement all the nonintronic regions
        //
        COMPLEMENT_NONINTRON (
            ch_not_intron_bed.map{
                meta, bed ->
                    meta.id = "intron"
                [ meta, bed ]},
            ch_chr_sizes_sorted)

        ch_intron_bed = COMPLEMENT_NONINTRON.out.bed
        ch_versions = ch_versions.mix(COMPLEMENT_NONINTRON.out.versions)

        //
        // MODULE: Convert the intron bed to gene pred
        //
        UCSC_BEDTOGENEPRED( ch_intron_bed )
        ch_genepred = UCSC_BEDTOGENEPRED.out.pred
        ch_versions = ch_versions.mix(UCSC_BEDTOGENEPRED.out.versions)

        //
        // MODULE: Convert the intron gene pred to gtf
        //
        UCSC_GENEPREDTOGTF( ch_genepred )
        ch_ucsc_gtf = UCSC_GENEPREDTOGTF.out.gtf
        ch_versions = ch_versions.mix(UCSC_GENEPREDTOGTF.out.versions)

        //
        // MODULE: Clean up the intron gtf
        //
        CREATE_INTRON_GTF ( ch_ucsc_gtf )
        ch_intron_gtf = CREATE_INTRON_GTF.out.intron_gtf
        ch_versions = ch_versions.mix(CREATE_INTRON_GTF.out.versions)

        // Cat the exon and intron gtf
        ch_gtfs = Channel.empty()
        ch_sorted_exon_gtf
            .map { meta, bed -> [bed] }
            .concat (
                ch_intron_gtf
                    .map { meta, bed -> [bed] }
                )
            .collect()
            .toList()
            .set{ ch_gtfs }

        //
        // MODULE: Combine all the gtfs together
        //
        ch_cat_gtf_in = Channel.of(["id": "exons_introns_merged"]).concat(ch_gtfs).collect()
        CAT_GTF (ch_cat_gtf_in)
        ch_versions = ch_versions.mix(CAT_GTF.out.versions)

        ch_prepared_gtf = Channel.empty()

        CAT_GTF.out.file_out
            .map {
                meta, gtf ->
                [ gtf ]
            }
            .set { ch_prepared_gtf }
        //ch_prepared_gtf = gtf

    } else {
        ch_prepared_gtf = gtf
    }

    emit:
    prepped_gtf = ch_prepared_gtf
    versions = ch_versions
}
