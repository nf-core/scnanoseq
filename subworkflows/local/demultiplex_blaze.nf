//
// Performs demultiplexing using BLAZE
//

// MODULES
include { BLAZE                                             } from '../../modules/local/blaze'
include { PREEXTRACT_FASTQ                                  } from '../../modules/local/preextract_fastq'
include { CORRECT_BARCODES                                  } from '../../modules/local/correct_barcodes'
include { SPLIT_FILE as SPLIT_FILE_BC_FASTQ                 } from "../../modules/local/split_file"
include { SPLIT_FILE as SPLIT_FILE_BC_CSV                   } from "../../modules/local/split_file"
include { CAT_CAT as CAT_CAT_PREEXTRACT                     } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_BARCODE                        } from "../../modules/nf-core/cat/cat/main"
include { PIGZ_COMPRESS                                     } from '../../modules/nf-core/pigz/compress/main'
include { PIGZ_UNCOMPRESS as PIGZ_UNCOMPRESS_BC             } from "../../modules/nf-core/pigz/uncompress/main"
include { PIGZ_UNCOMPRESS as PIGZ_UNCOMPRESS_FASTQ          } from "../../modules/nf-core/pigz/uncompress/main"

workflow DEMULTIPLEX_BLAZE {
    take:
        ch_trimmed_reads_combined           // channel: [ val(meta), path(trimmed_reads_combined) ]
        whitelist                            // channel: [ val(meta), path(whitelist) ]

    main:
        ch_versions = Channel.empty()
        ch_extracted_fastq = Channel.empty()
        ch_corrected_bc_info = Channel.empty()

        //
        // MODULE: Uncompress fastq.gz
        //

        //
        // MODULE: Generate whitelist
        //

        PIGZ_UNCOMPRESS_FASTQ ( ch_trimmed_reads_combined )

        ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS_FASTQ.out.versions)
        ch_trimmed_reads_combined_fastq = PIGZ_UNCOMPRESS_FASTQ.out.file

        //
        // MODULE: Unzip whitelist
        //

        // Unzip the whitelist if needed
        if (whitelist.extension == "gz"){

                PIGZ_UNCOMPRESS_BC ( [[:], whitelist] )

                ch_whitelist =
                        PIGZ_UNCOMPRESS_BC.out.file
                                .map {
                                        meta, whitelist ->
                                        [whitelist]
                                }

                ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS_BC.out.versions)
        } else {
                ch_whitelist = whitelist
        }

        BLAZE ( ch_trimmed_reads_combined_fastq, ch_whitelist )

        ch_putative_bc = BLAZE.out.putative_bc
        ch_gt_whitelist = BLAZE.out.whitelist
        ch_whitelist_bc_count = BLAZE.out.bc_count
        ch_versions = ch_versions.mix(BLAZE.out.versions)

        ch_split_bc_fastqs = ch_trimmed_reads_combined_fastq
        ch_split_bc = ch_putative_bc
        if (params.split_amount > 0) {
                SPLIT_FILE_BC_FASTQ( ch_trimmed_reads_combined_fastq, '.fastq', params.split_amount * 4 )

                SPLIT_FILE_BC_FASTQ.out.split_files
                        .transpose()
                        .set { ch_split_bc_fastqs }

                ch_versions = ch_versions.mix(SPLIT_FILE_BC_FASTQ.out.versions)

                SPLIT_FILE_BC_CSV ( ch_putative_bc, '.csv', (params.split_amount ) )
                SPLIT_FILE_BC_CSV.out.split_files
                        .transpose()
                        .set { ch_split_bc }
        }

        //
        // MODULE: Extract barcodes
        //

        PREEXTRACT_FASTQ( ch_split_bc_fastqs.join(ch_split_bc) )
        ch_barcode_info = PREEXTRACT_FASTQ.out.barcode_info
        ch_preextract_fastq = PREEXTRACT_FASTQ.out.extracted_fastq

        //
        // MODULE: Correct Barcodes
        //

        CORRECT_BARCODES (
                ch_barcode_info
                        .combine ( ch_gt_whitelist, by: 0)
                        .combine ( ch_whitelist_bc_count, by: 0 )
        )
        ch_corrected_bc_file = CORRECT_BARCODES.out.corrected_bc_info
        ch_versions = ch_versions.mix(CORRECT_BARCODES.out.versions)

        ch_extracted_fastq = ch_preextract_fastq
        ch_corrected_bc_info = ch_corrected_bc_file

        if (params.split_amount > 0){
                //
                // MODULE: Cat Preextract
                //
                CAT_CAT_PREEXTRACT(ch_preextract_fastq.groupTuple())
                ch_cat_preextract_fastq = CAT_CAT_PREEXTRACT.out.file_out

                //
                // MODULE: Cat barcode file
                //
                CAT_CAT_BARCODE (ch_corrected_bc_file.groupTuple())
                ch_corrected_bc_info = CAT_CAT_BARCODE.out.file_out

                //
                // MODULE: Zip the reads
                //
                PIGZ_COMPRESS (ch_cat_preextract_fastq )
                ch_extracted_fastq = PIGZ_COMPRESS.out.archive
                ch_versions = ch_versions.mix(PIGZ_COMPRESS.out.versions)
        }
    emit:
        // Versions
        versions = ch_versions

        extracted_fastq = ch_extracted_fastq
        corrected_bc_info = ch_corrected_bc_info
}
