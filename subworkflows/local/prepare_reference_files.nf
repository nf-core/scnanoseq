//
// Creates gtfs to that add introns as features
//

include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../../modules/nf-core/gunzip/main'
include { PREPARE_GTF            } from '../../subworkflows/local/prepare_gtf'

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

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
            ch_prepared_fasta = GUNZIP_FASTA( [ [:], fasta ]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = fasta  
        }

        ch_gtf = Channel.empty()
        if (gtf.endsWith('.gz')){
            ch_prepared_gtf = GUNZIP_GTF( [ [:], gtf ]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = gtf
        }
        
        //
        // MODULE: Index the fasta
        //
        SAMTOOLS_FAIDX( ch_prepared_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }

        //
        // SUBWORKFLOW: Prepare GTF
        //
        PREPARE_GTF (gtf_preparation_method, ch_gtf, fasta)
        ch_prepared_gtf = PREPARE_GTF.out.prepped_gtf
        ch_versions = ch_versions.mix(PREPARE_GTF.out.versions)


    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        prepped_gtf = ch_prepared_gtf
        versions = ch_versions
}
