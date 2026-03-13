process SEURAT {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-base conda-forge::r-seurat=4.1.1 conda-forge::r-ggplot2 conda-forge::r-optparse conda-forge::r-stringi"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b4cd78f1471a75cb2d338d3be506b2352723c0d2:4d30026c33d2809f4bf7b3a62f0e2b8529cb6915-0' :
        'biocontainers/mulled-v2-b4cd78f1471a75cb2d338d3be506b2352723c0d2:4d30026c33d2809f4bf7b3a62f0e2b8529cb6915-0' }"

    input:
    tuple val(meta), path(counts), path(flagstat)
    val mtx_format

    output:
    tuple val(meta), path("*.csv"), emit: seurat_stats
    tuple val(meta), path("*.png"), emit: seurat_qcs
    tuple val(meta), path("*.rds"), emit: seurat_rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (mtx_format.equals("MEX")) {
        """
        mkdir indir
        for file in $counts
        do
            mv \$file indir
        done

        seurat_qc.R \\
            $args \\
            -j indir \\
            -s $flagstat \\
            -d $prefix \\
            -r $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        END_VERSIONS
        """
    } else {
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
            r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        END_VERSIONS
        """
    }
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats.csv
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
