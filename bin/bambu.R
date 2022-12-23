#!/usr/bin/env Rscript

######################
### LOAD LIBRARIES ###
######################
library(optparse)
library(bambu)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
  make_option(c("-i", "--input_bam_dir" ), type="character", default=NULL    , metavar="path"   , help="Path to input bam file(s)."),
  make_option(c("-g", "--genome"        ), type="character", default=NULL    , metavar="path"   , help="Path to a fasta genome file."),
  make_option(c("-a", "--annotation"    ), type="character", default=NULL    , metavar="path"   , help="Path to gtf annotation file."),
  make_option(c("-o", "--outdir"        ), type="character", default="./"    , metavar="path"   , help="Output directory."),
  make_option(c("-r", "--outprefix"     ), type="character", default="bambu_", metavar="string" , help="Output prefix.")
)

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_bam_dir)) {
  print_help(opt_parser)
  stop("Please provide path to bam file(s).", call. = FALSE)
}

if (is.null(opt$genome)) {
  print_help(opt_parser)
  stop("Please provide the genome fasta file.", call. = FALSE)
}

if (is.null(opt$annotation)) {
  print_help(opt_parser)
  stop("Please provide the gtf annotation file file.", call. = FALSE)
}

################
### BAM PATH ###
################

# vector of path(s) to bam(s)
bam_files <- list.files(path = opt$input_bam_dir,
                        pattern = "\\.bam$",
                        full.names = TRUE)

bam_files

##########################
### PREPARE ANNOTATION ###
##########################

bambuAnnotations <- prepareAnnotations(opt$annotation)

############################
### BAMBU QUANTIFICATION ###
############################

#NOTE: trackReads = TRUE to map transcripts to reads
se <- bambu(reads = bam_files,
            annotations = bambuAnnotations,
            genome = opt$genome,
            trackReads = TRUE)

se

####################
### SAVE OUTPUTS ###
####################

# setting work dir to output dir
if (file.exists(opt$outdir) == FALSE) {
  dir.create(opt$outdir, recursive=TRUE)
}

setwd(opt$outdir)

# save standard bambu outputs
writeBambuOutput(se, path = "./", prefix = opt$outprefix)

# save bambu object to save trackReads info
saveRDS(se, file = "./scnanoseq_bambu.rds")

####################
### Session Info ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
