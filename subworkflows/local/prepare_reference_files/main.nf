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
        val_genome_fasta     // file: path/to/genome.fasta
        val_transcript_fasta // file: path/to/transcript.fasta
        val_gtf              // file: path/to/genome.gtf

    main:

        // Check if fasta and gtf are zipped

        //
        // MODULE: Unzip Genome FASTA
        //
        ch_prepared_genome_fasta = channel.empty()
        ch_prepared_genome_fai = channel.empty()

        if (val_genome_fasta) {
            if (val_genome_fasta.endsWith('.gz')){
                GUNZIP_GENOME_FASTA( [ [:], val_genome_fasta ])

                ch_prepared_genome_fasta = GUNZIP_GENOME_FASTA.out.file

            } else if (val_genome_fasta.endsWith('.zip')){
                UNZIP_GENOME_FASTA( [ [:], val_genome_fasta ])

                ch_prepared_genome_fasta = UNZIP_GENOME_FASTA.out.files

            } else {
                ch_prepared_genome_fasta = channel.of([ [:], val_genome_fasta ])
            }

            //
            // MODULE: Index the genome fasta
            //
            GENOME_FAIDX(
                ch_prepared_genome_fasta
                    .map {
                        meta, fasta ->
                        [meta, fasta, "$projectDir/assets/dummy_file.txt"]
                    },
                false
            )
            ch_prepared_genome_fai = GENOME_FAIDX.out.fai
        }

        //
        // MODULE: Unzip Transcript FASTA
        //
        ch_prepared_transcript_fasta = channel.empty()
        ch_prepared_transcript_fai = channel.empty()

        if (val_transcript_fasta) {
            if (val_transcript_fasta.endsWith('.gz')){
                GUNZIP_TRANSCRIPT_FASTA( [ [:], val_transcript_fasta ])

                ch_prepared_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.file

            } else if (val_transcript_fasta.endsWith('.zip')) {
                UNZIP_TRANSCRIPT_FASTA( [ [:], val_transcript_fasta ])

                ch_prepared_transcript_fasta = UNZIP_TRANSCRIPT_FASTA.out.files

            } else {
                ch_prepared_transcript_fasta = channel.of([ [:], val_transcript_fasta ])
            }

            //
            // MODULE: Index the transcript fasta
            //
            TRANSCRIPT_FAIDX(
                ch_prepared_transcript_fasta
                    .map {
                        meta, fasta ->
                        [meta, fasta, "$projectDir/assets/dummy_file.txt"]
                    },
                false
            )
            ch_prepared_transcript_fai = TRANSCRIPT_FAIDX.out.fai
        }

        //
        // MODULE: Unzip GTF
        //
        ch_meta_gtf = [ ['id': 'gtf'], val_gtf]
        ch_prepared_gtf = channel.empty()

        if (val_gtf.endsWith('.gz')){
            GUNZIP_GTF( ch_meta_gtf )

            ch_prepared_gtf = GUNZIP_GTF.out.file

        } else if (val_gtf.endsWith('.zip')) {
            UNZIP_GTF( ch_meta_gtf )

            ch_prepared_gtf = UNZIP_GTF.out.files

        } else {
            ch_prepared_gtf = ch_meta_gtf
        }

    emit:
        prepped_genome_fasta     = ch_prepared_genome_fasta     // channel: [ val(meta), path(genome_fasta) ]
        genome_fai               = ch_prepared_genome_fai       // channel: [ val(meta), path(genome_fai) ]
        prepped_transcript_fasta = ch_prepared_transcript_fasta // channel: [ val(meta), path(transcript_fasta) ]
        transcript_fai           = ch_prepared_transcript_fai   // channel: [ val(meta), path(transcript_fai) ]
        prepped_gtf              = ch_prepared_gtf              // channel: [ val(meta), path(gtf) ]
}
