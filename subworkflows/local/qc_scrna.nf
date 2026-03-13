//
// Performs feature quantification for long read single-cell rna data
//

include { SEURAT               } from '../../modules/local/seurat'
include { COMBINE_SEURAT_STATS } from '../../modules/local/combine_seurat_stats'

workflow QC_SCRNA {
    take:
        in_mtx
        in_flagstat
        mtx_format

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Seurat
        //
        SEURAT ( in_mtx.join(in_flagstat, by: [0]), mtx_format )
        ch_versions = ch_versions.mix(SEURAT.out.versions)

        //
        // MODULE: Combine Seurat Stats
        //
        COMBINE_SEURAT_STATS ( SEURAT.out.seurat_stats.collect{it[1]} )
        ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS.out.versions)

    emit:
        seurat_stats = COMBINE_SEURAT_STATS.out.combined_stats
        versions = ch_versions
}
