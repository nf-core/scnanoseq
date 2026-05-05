//
// Modifies the reference files for easier analysis
//

include { PIGZ_UNCOMPRESS as GUNZIP_GENOME_FASTA     } from '../../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as GUNZIP_GTF              } from '../../../modules/nf-core/pigz/uncompress/main'
include { UNZIPFILES as UNZIP_GENOME_FASTA           } from '../../../modules/nf-core/unzipfiles/main'
include { UNZIPFILES as UNZIP_TRANSCRIPT_FASTA       } from '../../../modules/nf-core/unzipfiles/main'
include { UNZIPFILES as UNZIP_GTF                    } from '../../../modules/nf-core/unzipfiles/main'

include { SAMTOOLS_FAIDX as GENOME_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as TRANSCRIPT_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        genome_fasta     // file: path/to/genome.fasta
        transcript_fasta // file: path/to/transcript.fasta
        gtf              // file: path/to/genome.gtf

    main:

        // Check if fasta and gtf are zipped

        //
        // MODULE: Unzip Genome FASTA
        //
        ch_genome_fasta = channel.empty()
        ch_genome_fai = channel.empty()

        if (genome_fasta) {
            if (genome_fasta.endsWith('.gz')){
                GUNZIP_GENOME_FASTA( [ [:], genome_fasta ])

                ch_genome_fasta = GUNZIP_GENOME_FASTA.out.file

            } else if (genome_fasta.endsWith('.zip')){
                UNZIP_GENOME_FASTA( [ [:], genome_fasta ])

                ch_genome_fasta = UNZIP_GENOME_FASTA.out.files

            } else {
                ch_genome_fasta = channel.of([ [:], genome_fasta ])
            }

            //
            // MODULE: Index the genome fasta
            //
            GENOME_FAIDX(
                ch_genome_fasta
                    .map {
                        meta, fasta ->
                        [meta, fasta, "$projectDir/assets/dummy_file.txt"]
                    },
                false
            )
            ch_genome_fai = GENOME_FAIDX.out.fai
        }

        //
        // MODULE: Unzip Transcript FASTA
        //
        ch_transcript_fasta = channel.empty()
        ch_transcript_fai = channel.empty()

        if (transcript_fasta) {
            if (transcript_fasta.endsWith('.gz')){
                GUNZIP_TRANSCRIPT_FASTA( [ [:], transcript_fasta ])

                ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.file

            } else if (transcript_fasta.endsWith('.zip')) {
                UNZIP_TRANSCRIPT_FASTA( [ [:], transcript_fasta ])

                ch_transcript_fasta = UNZIP_TRANSCRIPT_FASTA.out.files

            } else {
                ch_transcript_fasta = channel.of([ [:], transcript_fasta ])
            }

            //
            // MODULE: Index the transcript fasta
            //
            TRANSCRIPT_FAIDX(
                ch_transcript_fasta
                    .map {
                        meta, fasta ->
                        [meta, fasta, "$projectDir/assets/dummy_file.txt"]
                    },
                false
            )
            ch_transcript_fai = TRANSCRIPT_FAIDX.out.fai
        }

        //
        // MODULE: Unzip GTF
        //
        ch_prepared_gtf = channel.empty()
        if (gtf.endsWith('.gz')){
            GUNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = GUNZIP_GTF.out.file

        } else if (gtf.endsWith('.zip')) {
            UNZIP_GTF( [ [:], gtf ])

            ch_prepared_gtf = UNZIP_GTF.out.files

        } else {
            ch_prepared_gtf = [ [:], gtf]
        }

    emit:
        prepped_genome_fasta     = ch_genome_fasta
        genome_fai               = ch_genome_fai
        prepped_transcript_fasta = ch_transcript_fasta
        transcript_fai           = ch_transcript_fai
        prepped_gtf              = ch_prepared_gtf
}
