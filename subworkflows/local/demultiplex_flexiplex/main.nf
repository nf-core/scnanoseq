//
// Demultiplex reads using flexiplex
//

include { PIGZ_UNCOMPRESS                           } from '../../../modules/nf-core/pigz/uncompress/main'
include { SPLIT_SEQ                                 } from '../../../modules/local/split_seq'
include { FLEXIPLEX_DISCOVERY                       } from '../../../modules/local/flexiplex/discovery/main'
include { FLEXIPLEX_FILTER                          } from '../../../modules/local/flexiplex/filter/main'
include { FLEXIPLEX_ASSIGN                          } from '../../../modules/local/flexiplex/assign/main'
include { CAT_FASTQ                                 } from '../../../modules/nf-core/cat/fastq/main'
include { MERGEBARCODECOUNTS as MERGE_BARCODES      } from '../../../modules/local/mergebarcodecounts/main'

workflow DEMULTIPLEX_FLEXIPLEX {
    take:
        reads
        whitelist

    main:
        ch_versions = channel.empty()
        ch_flexiplex_fastq = channel.empty()
        ch_flexiplex_barcodes = channel.empty()

        // Unzip the whitelist if needed
        if (whitelist.extension == "gz"){

            PIGZ_UNCOMPRESS ( [[:], whitelist] )

            ch_whitelist =
                PIGZ_UNCOMPRESS.out.file
                    .map {
                        _meta, wl ->
                        [wl]
                    }

            ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions)
        } else {
            ch_whitelist = whitelist
        }


        flexiplex_input = reads
        if (params.split_amount > 0) {
            //
            // MODULE: Split reads into parts
            //
            SPLIT_SEQ (
                reads,
                '.fastq.gz',
                params.split_amount
            )

            ch_versions = ch_versions.mix(SPLIT_SEQ.out.versions_split_seq)

            // Transpose channel and add part to metadata
            flexiplex_input = SPLIT_SEQ.out.split_files
                .map { meta, rds ->
                    def newmeta = [splitcount: rds.size()]
                    [meta + newmeta, rds] }
                .transpose()
                .map { meta , rds ->
                    def part = (rds =~ /.*part_(\d+)\.fastq(?:\.gz)?$/)[0][1]
                    def newmeta = [part: part]
                    [meta + newmeta, rds] }
        }

        //
        // MODULE: Run flexiplex
        //

        FLEXIPLEX_DISCOVERY (
            flexiplex_input
        )

        ch_versions = ch_versions.mix(FLEXIPLEX_DISCOVERY.out.versions_flexiplex_discovery)


        //
        // Merge barcode counts if split
        //
        ch_barcodes = FLEXIPLEX_DISCOVERY.out.barcode_counts
        if (params.split_amount > 0 ) {

            ch_flexiplex_barcodes = FLEXIPLEX_DISCOVERY.out.barcode_counts
                .map { meta, barcode_counts ->
                    def key = groupKey(meta.subMap('id', 'single_end', 'cell_counts', 'type'), meta.splitcount)
                    [key, barcode_counts] }
                .groupTuple()

            MERGE_BARCODES (
                ch_flexiplex_barcodes
            )

            ch_versions = ch_versions.mix(MERGE_BARCODES.out.versions_mergebarcodecounts)
            ch_barcodes = MERGE_BARCODES.out.barcode_counts
        }

        //
        // MODULE: Filter flexiplex
        //

        FLEXIPLEX_FILTER (
            ch_barcodes,
            ch_whitelist
        )

        ch_versions = ch_versions.mix(FLEXIPLEX_FILTER.out.versions_flexiplex_filter)
        ch_corrected_bc_info = FLEXIPLEX_FILTER.out.barcodes

        // Merge the reads and barcodes channels
        flexiplex_input_barcodes = flexiplex_input
            .combine(ch_corrected_bc_info)
            .map { meta, rds, meta2, barcodes -> {
                meta.id == meta2.id ? [meta, rds, barcodes] : null }
            }

        //
        // MODULE: Assign flexiplex
        //

        FLEXIPLEX_ASSIGN (
            flexiplex_input_barcodes,
        )

        ch_versions = ch_versions.mix(FLEXIPLEX_ASSIGN.out.versions_flexiplex_assign)

        ch_flexiplex_fastq = FLEXIPLEX_ASSIGN.out.reads
        if (params.split_amount > 0) {
            //
            // MODULE: cat fastq
            //

            ch_grouped_flexiplex_fastq = FLEXIPLEX_ASSIGN.out.reads
            .map { meta, rds ->
                def key = groupKey(meta.subMap('id', 'single_end', 'cell_counts', 'type'), meta.splitcount)
                [key, rds] }
            .groupTuple()

            CAT_FASTQ (
                ch_grouped_flexiplex_fastq
            )

            ch_flexiplex_fastq = CAT_FASTQ.out.reads
        }

    emit:
        flexiplex_fastq = ch_flexiplex_fastq
        flexiplex_barcodes = ch_corrected_bc_info
        versions = ch_versions
}
