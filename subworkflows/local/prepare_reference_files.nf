//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNZIP_GTF   } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta_preparation_method
        gtf_preparation_method
        fasta
        gtf

    main:
        ch_versions = Channel.empty()

        // Check if fasta and gtf are zipped

        //
        ch_prepared_fasta = Channel.empty()
        if (fasta.endsWith('.gz')){
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.gunzip
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = [ [:], fasta ]
        }

        ch_prepared_gtf = Channel.empty()
        if (gtf.endsWith('.gz')){
            UNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = UNZIP_GTF.out.gunzip
            ch_versions = ch_versions.mix(UNZIP_GTF.out.versions)
        } else {
            ch_prepared_gtf = [ [:], gtf]
        }

        //
        // MODULE: Index the fasta
        //
        SAMTOOLS_FAIDX( ch_prepared_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        prepped_gtf = ch_prepared_gtf
        versions = ch_versions
}
