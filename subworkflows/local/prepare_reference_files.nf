//
// Creates gtfs to that add introns as features
//

include { PREPARE_GTF } from '../../subworkflows/local/prepare_gtf'

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

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

        //
        // MODULE: Index the fasta
        //
        SAMTOOLS_FAIDX( [ [ "id": "fasta"], fasta ], [ ["id": "fai"], "$projectDir/assets/dummy_file.txt" ])
        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        prepped_gtf = ch_prepared_gtf
        versions = ch_versions
}
