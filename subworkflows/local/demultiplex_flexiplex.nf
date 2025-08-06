//
// Creates gtfs to that add introns as features
//

include { PIGZ_COMPRESS                             } from '../../modules/nf-core/pigz/compress/main'
include { SEQKIT_SPLIT2                             } from '../../modules/nf-core/seqkit/split2/main'
include { FLEXIPLEX_DISCOVERY                       } from '../../modules/local/flexiplex/discovery/main'
include { FLEXIPLEX_FILTER                          } from '../../modules/local/flexiplex/filter/main'
include { FLEXIPLEX_ASSIGN                          } from '../../modules/local/flexiplex/assign/main'
include { CAT_FASTQ                                 } from '../../modules/nf-core/cat/fastq/main'


workflow DEMULTIPLEX_FLEXIPLEX {
    take:
        ch_reads
        whitelist

    main:
        ch_versions = Channel.empty()
        ch_flexiplex_fastq = Channel.empty()
        ch_flexiplex_barcodes = Channel.empty()
        
        //
        // MODULE: Run flexiplex
        //
        FLEXIPLEX_DISCOVERY (
            ch_reads
    	)
        
        ch_versions = ch_versions.mix(FLEXIPLEX_DISCOVERY.out.versions)
        
        // 
        // MODULE: Filter flexiplex
        //
        
        FLEXIPLEX_FILTER (
            FLEXIPLEX_DISCOVERY.out.barcode_counts,
            whitelist
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_FILTER.out.versions)
        ch_corrected_bc_info = FLEXIPLEX_FILTER.out.barcodes
        
        //
        // MODULE: Run SEQKIT_SPLIT2
        //
        SEQKIT_SPLIT2 (
            ch_reads
        )
        
        // Transpose channel and add part to metadata
        SEQKIT_SPLIT2.out.reads
            | transpose
            | map { meta, reads ->
                part = (reads =~ /.*part_(\d+)\.fastq(?:\.gz)?$/)[0][1]
                newmap = [part: part]
                [meta + newmap, reads] }
            | set { ch_split_fastq }
                
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
        
        // Merge the reads and barcodes channels
        ch_split_fastq 
            | combine(ch_corrected_bc_info)
            | map { meta, reads, meta2, barcodes -> { 
                meta.id == meta2.id ? [meta, reads, barcodes] : null }}
            | set { ch_split_fastq_barcode }
    
        
        // 
        // MODULE: Assign flexiplex
        //
        FLEXIPLEX_ASSIGN (
            ch_split_fastq_barcode,
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_ASSIGN.out.versions)
        
        //
        // MODULE: Compress Fastqs
        //
        PIGZ_COMPRESS ( FLEXIPLEX_ASSIGN.out.reads )
        
        ch_versions = ch_versions.mix(PIGZ_COMPRESS.out.versions)
        
        // Group by ID for CATFASTQ
        PIGZ_COMPRESS.out.archive
            | map { meta, reads ->
                [meta.subMap('id', 'single_end'), meta.part, reads] }
            | groupTuple
            | map { meta, part, reads -> [meta + [partcount: part.size()], reads] }
            | set { ch_grouped_flexiplex_fastq }
        
        
        //
        // MODULE: cat fastq
        //
        
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