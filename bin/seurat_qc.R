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
    make_option(c("-d", "--id"            ), type="character", default="scnanoseq", metavar="integer", help="Project name for Seurat object."                                     ),
    make_option(c("-o", "--outdir"        ), type="character", default="./"       , metavar="path"   , help="Output directory."                                                   ),
    make_option(c("-r", "--outprefix"     ), type="character", default="seurat_qc", metavar="string" , help="Output prefix."                                                      )
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

######################
### READ IN INPUTs ###
######################

# cell or nuclei matrix (calling it cell for simplicity)
cell_bc_matrix <- read.table(opt$input_matrix, sep="\t", header = TRUE, row.names = 1)

flagstat_lines <- readLines(opt$flagstat)

#######################
### GET TOTAL READS ###
#######################

# This will grep the lines to get the line index containing the total read count
index_nums <- grep("total", flagstat_lines)

# This will parse out the total read count
total_reads <- as.numeric(gsub("([0-9]+).*$", "\\1", flagstat_lines[index_nums]))

#####################
### SEURAT OBJECT ###
#####################

# Create the Seurat object
#NOTE: we do not perform any pre-filtering at this point
seurat_obj <- CreateSeuratObject(counts = cell_bc_matrix,
                                 min.cells = 0,
                                 min.features = 0,
                                 project = opt$id)

######################
### GENERATE PLOTS ###
######################

# This is a set of very preliminary QCs from the input matrices

# Setting work dir to output dir
if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive=TRUE)
}

setwd(opt$outdir)

# Use prefix for out files
out_sample_qc_figs <- paste0(opt$outprefix,".png")

#TODO: again reminder to add density plot into this mix

# Generate the violin plot
violing_plot <- VlnPlot(object = seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), assay = "RNA")

# Feature Scatter plot
feature_plot <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

png(out_sample_qc_figs, width = 800, height = 800)

plot(violing_plot + feature_plot)

dev.off()


######################
### GENERATE STATS ###
######################

# Main goal here is to provide detailed information to users in the
# MultiQC report linked to read depth per cell, # of genes etc.

# Get the "Estimated Number of Cells"
est_cell_number <- length(colnames(cell_bc_matrix))

# Get the "Mean Reads per Cell"
mean_reads_per_cell <- round(total_reads / est_cell_number, digits = 2)

# Get the "Median Genes per Cell"
median_genes_per_cell <- median(seurat_obj$nFeature_RNA)

####################
### OUTPUT TABLE ###
####################

# Acquire all stats
output_table <- data.frame(est_cell_number, mean_reads_per_cell, median_genes_per_cell)

# Set the colnames for the output table
colnames(output_table) <- c("Estimated Cell Number", "Mean Reads Per Cell", "Median Genes Per Cell")

# Write out the data frame
out_stats <- paste0(opt$outprefix,".csv")

write.csv(output_table, out_stats, row.names = FALSE)

####################
### Session Info ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
