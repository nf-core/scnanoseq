#!/usr/bin/env Rscript

######################
### LOAD LIBRARIES ###
######################
library(optparse)
library(Seurat)

# TODO: add density plot custom function for more QC

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(c("-i", "--input_matrix"  ), type="character", default=NULL       , metavar="path"   , help="Count file matrix where rows are genes and columns are cells/nuclei."),
    make_option(c("-s", "--flagstat"      ), type="character", default=NULL       , metavar="path"   , help="Flagstat file from samtools QC."                                     ),
    make_option(c("-c", "--min_cells"     ), type="integer"  , default=0          , metavar="integer", help="Minimun number of cells to retain during CreateSeuratObject."        ),
    make_option(c("-f", "--min_features"  ), type="integer"  , default=0          , metavar="integer", help="Minimun number of features to retain during CreateSeuratObject."     ),
    make_option(c("-d", "--id"            ), type="character", default="scnanoseq", metavar="integer", help="Project name for Seurat object."                                     ),
    make_option(c("-o", "--outdir"        ), type="character", default="./"       , metavar="path"   , help="Output directory."                                                   ),
    make_option(c("-p", "--outprefix"     ), type="character", default="seurat_qc", metavar="string" , help="Output prefix."                                                      )
)

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_matrix)) {
    print_help(opt_parser)
    stop("Please provide a single-cell/nuclei matrix.", call. = FALSE)
}

if (is.null(opt$flagstat)) {
    print_help(opt_parser)
    stop("Please provide the samtools flagstat file.", call. = FALSE)
}
