//
// Creates gtfs to that add introns as features
//

include { TRANSCRIPT_TO_EXON } from '../../modules/local/prepare_gtf_transcript_to_exon'

workflow PREPARE_REFERENCE_FILES {
    take:
    gtf_preparation_method
    gtf

    main:

    if (gtf_preparation_method == "1") {
        // Convert the Transcript to exons
        TRANSCRIPT_TO_EXON(gtf)
        ch_prepared_gtf = TRANSCRIPT_TO_EXON.out.ch_processed_gtf

    } else if (gtf_preparation_method == "2") {
        // Get the chromosome sizes
        
        // Sort the gtf file
        
        // Complement the gtf to the chr sizes to get the intergenic regions
        
        // Get the exon regions
        
        // Get the intron regions
        
        // Generate the intron gtf
        
        // Clean up the intron gtf
        
        // Cat the the exon and intron gtf
    } else {
        ch_prepared_gtf = gtf
    }

    emit:
    ch_prepared_gtf
}
