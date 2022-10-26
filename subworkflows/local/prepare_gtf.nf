//
// Creates gtfs to that add introns as features
//

include { TRANSCRIPT_TO_EXON                          } from '../../modules/local/prepare_gtf_transcript_to_exon'
include { SORT_GTF                                    } from '../../modules/local/sort_gtf'
include { SORT_GTF as SORT_EXON_GTF                   } from '../../modules/local/sort_gtf'
include { GET_GTF_FEATURES                            } from '../../modules/local/get_gtf_features'
include { GTF2BED                                     } from '../../modules/local/gtf2bed'
include { UCSC_BEDTOGENEPRED                          } from '../../modules/local/ucsc_bedtogenepred'
include { UCSC_GENEPREDTOGTF                          } from '../../modules/local/ucsc_genepredtogtf'

include { CUSTOM_GETCHROMSIZES                        } from '../../modules/nf-core/custom/getchromsizes/main'
include { BEDTOOLS_COMPLEMENT as COMPLEMENT_GTF       } from '../../modules/nf-core/bedtools/complement/main'
include { BEDTOOLS_COMPLEMENT as COMPLEMENT_NONINTRON } from '../../modules/nf-core/bedtools/complement/main'
include { CAT_CAT                                     } from '../../modules/nf-core/cat/cat/main'

workflow PREPARE_GTF {
    take:
    gtf_preparation_method
    gtf
    fasta

    main:

    if (gtf_preparation_method == "1") {
        // Convert the Transcript to exons
        TRANSCRIPT_TO_EXON(gtf)
        ch_prepared_gtf = TRANSCRIPT_TO_EXON.out.ch_processed_gtf

    } else if (gtf_preparation_method == "2") {
        ch_prepared_gtf = gtf
        // Get the chromosome sizes
        CUSTOM_GETCHROMSIZES ( [ ["id":"chr_sizes"], fasta ])
        chr_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { meta, chr_sizes -> [chr_sizes] }

        // Sort the gtf
        SORT_GTF ( [ [ "id": "base" ], gtf ])
        ch_sorted_gtf = SORT_GTF.out.gtf
        
        // Complement the gtf to the chr sizes to get the intergenic regions
        ch_intergenic_in = Channel.empty()
        ch_sorted_gtf
            .map {
                meta, gtf ->
                    meta.id = "intergenic"
                    [ meta, gtf ]
            }
            .set{ch_intergenic_in}

        COMPLEMENT_GTF ( ch_intergenic_in,
                              chr_sizes )
        ch_intergenic_bed = COMPLEMENT_GTF.out.bed
        
        // Get the exon regions
        GET_GTF_FEATURES( [ [ "id": "exon" ], gtf], "exons" )
        ch_exon_gtf = GET_GTF_FEATURES.out.gtf

        SORT_EXON_GTF(ch_exon_gtf)
        ch_sorted_exon_gtf = SORT_EXON_GTF.out.gtf

        GTF2BED( ch_sorted_exon_gtf )
        ch_exon_bed = GTF2BED.out.bed


        // Get the intron regions
        
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

        // Now we combine the file list with a meta object
        ch_cat_files_in = Channel.of(["id": "not_introns"]).concat(ch_files).collect()
        
        CAT_CAT( ch_cat_files_in)

        ch_not_intron_bed = CAT_CAT.out.file_out

        COMPLEMENT_NONINTRON ( 
            ch_not_intron_bed
                .map { meta, gtf -> 
                    meta.id = "intron"
                    [ meta, gtf ]
                },
            chr_sizes)

        ch_intron_bed = COMPLEMENT_NONINTRON.out.bed

        // Generate the intron gtf

        
        // Clean up the intron gtf
        
        // Cat the the exon and intron gtf
    } else {
        ch_prepared_gtf = gtf
    }

    emit:
    ch_prepared_gtf
}
