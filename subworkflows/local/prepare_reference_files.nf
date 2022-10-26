//
// Creates gtfs to that add introns as features
//

include { PREPARE_GTF } from '../../subworkflows/local/prepare_gtf'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta_preparation_method
        gtf_preparation_method
        fasta
        gtf

    main:
        // There's nothing we wish to change about the fasta, so just need to return it
        ch_prepared_fasta = fasta
        
        //
        // SUBWORKFLOW: Prepare GTF
        //
        PREPARE_GTF (gtf_preparation_method, gtf, fasta)
        ch_prepared_gtf = PREPARE_GTF.out.ch_prepared_gtf

    emit:
        ch_prepared_fasta
        ch_prepared_gtf
}
