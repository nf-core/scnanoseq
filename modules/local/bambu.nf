process BAMBU {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-bambu=3.0.5 conda-forge::r-optparse" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-451c38e14eb84db642390085de685f74f02238cb:1c7bc0278b2279671d4eed26b47a68576a78b0a7-0' :
        'quay.io/biocontainers/mulled-v2-451c38e14eb84db642390085de685f74f02238cb:1c7bc0278b2279671d4eed26b47a68576a78b0a7-0' }"

    input:
    tuple path(fasta), path(gtf)
    path bams

    output:
    path "*counts_transcript.txt"    , emit: ch_transcript_counts
    path "*extended_annotations.gtf" , emit: extended_gtf
    path "*.rds"                     , emit: rds
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bams = bams.join(",")

    //TODO add cores option below once implemented in the script

    """
    bambu.R \\
        $args \\
        -i $bams \\
        -g $fasta \\
        -a $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
    END_VERSIONS
    """
}
