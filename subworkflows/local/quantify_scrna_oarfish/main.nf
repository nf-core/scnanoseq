//
// Performs feature quantification for long read single-cell rna data
//

include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { OARFISH       } from '../../../modules/local/oarfish'
include { QC_SCRNA      } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_OARFISH {
    take:
        ch_bam      // channel: [ val(meta), path(bam) ]
        ch_flagstat // channel: [ val(meta), path(flagstat) ]
        ch_fasta    // channel: [ val(meta), path(fasta) ]

        val_skip_qc     // bool: Skip qc steps
        val_skip_seurat // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        //
        // MODULE: Samtools Sort
        //
        SAMTOOLS_SORT ( ch_bam, ch_fasta.first(), "bai" )

        //
        // MODULE: Oarfish
        //
        OARFISH ( SAMTOOLS_SORT.out.bam )
        ch_versions = ch_versions.mix(OARFISH.out.versions_oarfish)

        ch_transcript_qc_stats = channel.empty()
        if (!val_skip_qc && !val_skip_seurat) {
            QC_SCRNA(
                OARFISH.out.features
                    .join(OARFISH.out.barcodes, by: [0])
                    .join(OARFISH.out.mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            def new_meta = [ 'id' : meta.id ]
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA.out.versions)
        }

    emit:
        versions            = ch_versions            // channel: [ val(meta), path(versions) ]
        transcript_qc_stats = ch_transcript_qc_stats // channel: [ val(meta), path(transcript_qc_stats) ]
}
