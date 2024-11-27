//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA         } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNZIP_GTF           } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                         } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SPLIT } from '../../modules/nf-core/samtools/faidx/main'
include { SPLIT_GTF                              } from '../../modules/local/split_gtf'
include { SPLIT_FASTA                            } from '../../modules/local/split_fasta'

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

            ch_prepared_fasta = UNZIP_FASTA.out.file
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = [ [:], fasta ]
        }

        ch_prepared_gtf = Channel.empty()
        if (gtf.endsWith('.gz')){
            UNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = UNZIP_GTF.out.file
            ch_versions = ch_versions.mix(UNZIP_GTF.out.versions)
        } else {
            ch_prepared_gtf = [ [:], gtf]
        }

        //
        // MODULE: Index the fasta
        //
        SAMTOOLS_FAIDX( ch_prepared_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai

        //
        // MODULE: Split the FASTA
        //
        SPLIT_FASTA( ch_prepared_fasta )
        ch_split_fasta = SPLIT_FASTA.out.split_fasta
            .flatten()
            .map{
                fasta ->
                    fasta_basename = fasta.toString().split('/')[-1]
                    meta = [ 'chr': fasta_basename.split(/\./)[0] ]
                    [ meta, fasta ]
            }

        SAMTOOLS_FAIDX_SPLIT( ch_split_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_split_fai = SAMTOOLS_FAIDX_SPLIT.out.fai

        //
        // MODULE: Split the GTF
        //
        SPLIT_GTF( ch_prepared_gtf )
        ch_split_gtf = SPLIT_GTF.out.split_gtf
            .flatten()
            .map{
                gtf ->
                    gtf_basename = gtf.toString().split('/')[-1]
                    meta = ['chr': gtf_basename.split(/\./)[0]]
                    [ meta, gtf ]
            }

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        prepped_gtf = ch_prepared_gtf
        split_gtf = ch_split_gtf
        split_fasta = ch_split_fasta
        split_fai = ch_split_fai
        versions = ch_versions
}
