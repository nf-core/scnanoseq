//
// Performs feature quantification for long read single-cell rna data
//

include { SEURAT               } from '../../../modules/local/seurat'
include { COMBINE_SEURAT_STATS } from '../../../modules/local/combine_seurat_stats'

workflow QC_SCRNA {
    take:
        ch_mtx         // channel: [ val(meta), path(cell_bc_matrix) ]
        ch_flagstat    // channel: [ val(meta), path(flagstat) ]
        val_mtx_format // str: Format of the cell barcode matrix (e.g. "MEX", "MTX")

    main:
        ch_versions = channel.empty()

        //
        // MODULE: Seurat
        //
        SEURAT ( ch_mtx.join(ch_flagstat, by: [0]), val_mtx_format )
        ch_versions = ch_versions.mix(SEURAT.out.versions_seurat)

        //
        // MODULE: Combine Seurat Stats
        //
        COMBINE_SEURAT_STATS ( SEURAT.out.seurat_stats.collect{ v -> v[1] } )
        ch_versions = ch_versions.mix(COMBINE_SEURAT_STATS.out.versions_combine_seurat_stats)

    emit:
        seurat_stats = COMBINE_SEURAT_STATS.out.combined_stats // channel: [ val(meta), path(seurat_stats) ]
        versions = ch_versions                                 // channel: [ val(meta), path(versions) ]
}
