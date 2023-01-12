process BAMBU {
    label 'process_medium'

    conda ("bioconda::bioconductor-bambu=3.0.5 conda-forge::r-optparse")
    // TODO: add mulled container once PR is merged (PR already present in biocontainers repo)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //    'quay.io/biocontainers/YOUR-TOOL-HERE' }"

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
