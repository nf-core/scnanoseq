process SEURAT {
    tag "$meta.id"
    label 'process_low'

    // TODO: figure out container for this one ; will note the following for now
    // seurat-scripts:4.0.0--hdfd78af_0 (need to test and see if it contains anything else needed)
    // or find alternatives
    conda ("bioconda::r-seurat=4.1.1")

    input:
    tuple val(meta), path(counts), path(flagstat)

    output:
    tuple val(meta), path("*.csv"), emit: seurat_stats
    tuple val(meta), path("*.png"), emit: seurat_qcs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seurat_qc.R \\
        $args \\
        -i $counts \\
        -s $flagstat \\
        -d $prefix \\
        -r $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
