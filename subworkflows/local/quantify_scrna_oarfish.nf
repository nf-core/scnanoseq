//
// Performs feature quantification for long read single-cell rna data
//

include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { OARFISH       } from '../../modules/local/oarfish'
include { QC_SCRNA      } from '../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_OARFISH {
    take:
        in_bam
        in_bai
        in_flagstat
        in_fasta
        skip_qc
        skip_seurat

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Samtools Sort
        //
        SAMTOOLS_SORT ( in_bam, in_fasta )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

        //
        // MODULE: Oarfish
        //
        OARFISH ( SAMTOOLS_SORT.out.bam )
        ch_versions = ch_versions.mix(OARFISH.out.versions)

        ch_transcript_qc_stats = Channel.empty()
        if (!params.skip_qc && !params.skip_seurat) {
            QC_SCRNA(
                OARFISH.out.features
                    .join(OARFISH.out.barcodes, by: [0])
                    .join(OARFISH.out.mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            new_meta = [ 'id' : meta.id ]
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                in_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA.out.versions)
        }

    emit:
        versions            = ch_versions
        transcript_qc_stats = ch_transcript_qc_stats
}
