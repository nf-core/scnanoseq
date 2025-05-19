//
// Modifies the reference files for easier analysis
//

include { PIGZ_UNCOMPRESS as GUNZIP_GENOME_FASTA     } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as GUNZIP_GTF              } from '../../modules/nf-core/pigz/uncompress/main'
include { UNZIPFILES as UNZIP_GENOME_FASTA           } from '../../modules/nf-core/unzipfiles/main'
include { UNZIPFILES as UNZIP_TRANSCRIPT_FASTA       } from '../../modules/nf-core/unzipfiles/main'
include { UNZIPFILES as UNZIP_GTF                    } from '../../modules/nf-core/unzipfiles/main'

include { SAMTOOLS_FAIDX as GENOME_FAIDX            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as TRANSCRIPT_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        genome_fasta     // file: path/to/genome.fasta
        transcript_fasta // file: path/to/transcript.fasta
        gtf              // file: path/to/genome.gtf

    main:
        ch_versions = Channel.empty()

        // Check if fasta and gtf are zipped

        //
        // MODULE: Unzip Genome FASTA
        //
        ch_genome_fasta = Channel.empty()
        ch_genome_fai = Channel.empty()

        if (genome_fasta) {
            if (genome_fasta.endsWith('.gz')){
                GUNZIP_GENOME_FASTA( [ [:], genome_fasta ])

                ch_genome_fasta = GUNZIP_GENOME_FASTA.out.file
                ch_versions = ch_versions.mix(GUNZIP_GENOME_FASTA.out.versions)

            } else if (genome_fasta.endsWith('.zip')){
                UNZIP_GENOME_FASTA( [ [:], genome_fasta ])

                ch_genome_fasta = UNZIP_GENOME_FASTA.out.files
                ch_versions = ch_versions.mix(UNZIP_GENOME_FASTA.out.versions)

            } else {
                ch_genome_fasta = [ [:], genome_fasta ]
            }

            //
            // MODULE: Index the genome fasta
            //
            GENOME_FAIDX( ch_genome_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
            ch_genome_fai = GENOME_FAIDX.out.fai
        }

        //
        // MODULE: Unzip Transcript FASTA
        //
        ch_transcript_fasta = Channel.empty()
        ch_transcript_fai = Channel.empty()
        if (transcript_fasta) {
            if (transcript_fasta.endsWith('.gz')){
                GUNZIP_TRANSCRIPT_FASTA( [ [:], transcript_fasta ])

                ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.file
                ch_versions = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)

            } else if (transcript_fasta.endsWith('.zip')) {
                UNZIP_TRANSCRIPT_FASTA( [ [:], transcript_fasta ])

                ch_transcript_fasta = UNZIP_TRANSCRIPT_FASTA.out.files
                ch_versions = ch_versions.mix(UNZIP_TRANSCRIPT_FASTA.out.versions)

            } else {
                ch_transcript_fasta = [ [:], transcript_fasta ]
            }

            //
            // MODULE: Index the transcript fasta
            //
            TRANSCRIPT_FAIDX( ch_transcript_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
            ch_transcript_fai = TRANSCRIPT_FAIDX.out.fai
        }

        //
        // MODULE: Unzip GTF
        //
        ch_prepared_gtf = Channel.empty()
        if (gtf.endsWith('.gz')){
            GUNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = GUNZIP_GTF.out.file
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        } else if (gtf.endsWith('.zip')) {
            UNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = UNZIP_GTF.out.files
            ch_versions = ch_versions.mix(UNZIP_GTF.out.versions)

        } else {
            ch_prepared_gtf = [ [:], gtf]
        }

    emit:
        prepped_genome_fasta     = ch_genome_fasta
        genome_fai               = ch_genome_fai
        prepped_transcript_fasta = ch_transcript_fasta
        transcript_fai           = ch_transcript_fai
        prepped_gtf              = ch_prepared_gtf
        versions                 = ch_versions
}
