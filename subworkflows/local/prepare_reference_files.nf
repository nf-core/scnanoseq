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
        ch_versions = Channel.empty()

        // There's nothing we wish to change about the fasta, so just need to return it
        ch_prepared_fasta = fasta

        //
        // SUBWORKFLOW: Prepare GTF
        //
        PREPARE_GTF (gtf_preparation_method, gtf, fasta)
        ch_prepared_gtf = PREPARE_GTF.out.prepped_gtf
        ch_versions = ch_versions.mix(PREPARE_GTF.out.versions)

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_gtf = ch_prepared_gtf
        versions = ch_versions
}
