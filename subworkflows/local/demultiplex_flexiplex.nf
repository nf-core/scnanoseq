//
// Creates gtfs to that add introns as features
//

include { PIGZ_COMPRESS                             } from '../../modules/nf-core/pigz/compress/main'
include { PIGZ_UNCOMPRESS                           } from '../../modules/nf-core/pigz/uncompress/main'
include { SEQKIT_SPLIT2                             } from '../../modules/nf-core/seqkit/split2/main'
include { FLEXIPLEX_DISCOVERY                       } from '../../modules/local/flexiplex/discovery/main'
include { FLEXIPLEX_FILTER                          } from '../../modules/local/flexiplex/filter/main'
include { FLEXIPLEX_ASSIGN                          } from '../../modules/local/flexiplex/assign/main'
include { CAT_FASTQ                                 } from '../../modules/nf-core/cat/fastq/main'
include { MERGEBARCODECOUNTS as MERGE_BARCODES      } from '../../modules/local/mergebarcodecounts/main'

workflow DEMULTIPLEX_FLEXIPLEX {
    take:
        reads
        whitelist

    main:
        ch_versions = Channel.empty()
        ch_flexiplex_fastq = Channel.empty()
        ch_flexiplex_barcodes = Channel.empty()
        
        // Unzip the whitelist if needed
        if (whitelist.extension == "gz"){

            PIGZ_UNCOMPRESS ( [[:], whitelist] )

            ch_whitelist =
                PIGZ_UNCOMPRESS.out.file
                    .map {
                        meta, whitelist ->
                        [whitelist]
                    }

            ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions)
        } else {
            ch_whitelist = whitelist
        }
        
        //
        // MODULE: SPLIT2: Split reads into parts
        //
        SEQKIT_SPLIT2 (
            reads
        )
        
        ch_versions = SEQKIT_SPLIT2.out.versions
        
        // Transpose channel and add part to metadata
        ch_split_fastq = SEQKIT_SPLIT2.out.reads
            .map { meta, reads -> 
                  newmeta = [splitcount: reads.size()]
                  [meta + newmeta, reads] }
            .transpose()
            .map { meta , reads -> 
                  part = (reads =~ /.*part_(\d+)\.fastq(?:\.gz)?$/)[0][1]
                  newmeta = [part: part]
                  [meta + newmeta, reads] }
      
        //
        // MODULE: Run flexiplex
        //
        
        FLEXIPLEX_DISCOVERY (
            ch_split_fastq
    	)
        
        ch_versions = ch_versions.mix(FLEXIPLEX_DISCOVERY.out.versions)
        
        //
        // Merge barcode counts
        //
        
        ch_flexiplex_barcodes = FLEXIPLEX_DISCOVERY.out.barcode_counts
            .map { meta, barcode_counts -> 
                  key = groupKey(meta.subMap('id', 'single_end', 'cell_counts', 'type'), meta.splitcount)
                  [key, barcode_counts] }
            .groupTuple()

        MERGE_BARCODES (
            ch_flexiplex_barcodes
        )

        ch_versions = ch_versions.mix(MERGE_BARCODES.out.versions)

        //
        // MODULE: Filter flexiplex
        //
        
        FLEXIPLEX_FILTER (
            MERGE_BARCODES.out.barcode_counts,
            ch_whitelist
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_FILTER.out.versions)
        ch_corrected_bc_info = FLEXIPLEX_FILTER.out.barcodes
        
        // Merge the reads and barcodes channels
        ch_split_fastq_barcode = ch_split_fastq 
            .combine(ch_corrected_bc_info)
            .map { meta, reads, meta2, barcodes -> { 
                meta.id == meta2.id ? [meta, reads, barcodes] : null }}

        //
        // MODULE: Assign flexiplex
        //
        
        FLEXIPLEX_ASSIGN (
            ch_split_fastq_barcode,
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_ASSIGN.out.versions)
        
        //
        // MODULE: cat fastq
        //
        
        ch_grouped_flexiplex_fastq = FLEXIPLEX_ASSIGN.out.reads
            .map { meta, reads -> 
                  key = groupKey(meta.subMap('id', 'single_end', 'cell_counts', 'type'), meta.splitcount)
                  [key, reads] }
            .groupTuple()
        
        CAT_FASTQ (
            ch_grouped_flexiplex_fastq
        )
        
        ch_flexiplex_fastq = CAT_FASTQ.out.reads
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
        
    emit:
        flexiplex_fastq = ch_flexiplex_fastq
        flexiplex_barcodes = ch_corrected_bc_info
        versions = ch_versions
}