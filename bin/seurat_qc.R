#!/usr/bin/env Rscript

######################
### LOAD LIBRARIES ###
######################
library(optparse)
library(Seurat)
library(ggplot2)

#######################
### CUSTOM FUNCTION ###
#######################

#' plotSingleCellDensity
#' This function creates a density plot of the distribution of cells (nFeature) or molecules/UMIs (nCount)
#' @param input_obj A Seurat object
#' @param metadata_variable Metadata column name linked to cell numbers and UMIs - typically `nFeature_assay_name` or `nCount_assay_name`
#' @param group.by Metadata column name linked to the grouping variable to group the plot by
#' @param scale_x_axis Transforms x-axis to the log10 scale
#' @param geom_density_level Set transparency level for density plot - lower values generate more transparent density curves
#'
#' @return A density plot

plotSingleCellDensity <- function(input_obj,
                                    metadata_variable,
                                    group.by = "orig.ident",
                                    scale_x_axis = FALSE,
                                    geom_density_level = 0.2) {

    metadata <- dplyr::select(
        input_obj@meta.data,
        {{ metadata_variable }},
        {{ group.by }}
    )

    meta_density <- ggplot2::ggplot(
        metadata,
        aes(
            x = .data[[metadata_variable]],
            color = .data[[group.by]],
            fill = .data[[group.by]]
        )
    ) +
    geom_density(alpha = geom_density_level) +
    theme_classic()

    if (scale_x_axis == TRUE) {
        return(meta_density + scale_x_log10())
    } else {
        return(meta_density)
    }
}

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(c("-i", "--input_matrix"  ), type="character", default=NULL       , metavar="path"   , help="Count file matrix where rows are genes and columns are cells/nuclei."),
    make_option(c("-j", "--input_dir"     ), type="character", default=NULL       , metavar="path"   , help="Directory containing matrix.mtx, genes.tsv (or features.tsv) , and barcodes.tsv."),
    make_option(c("-s", "--flagstat"      ), type="character", default=NULL       , metavar="path"   , help="Flagstat file from samtools QC."                                     ),
    make_option(c("-d", "--id"            ), type="character", default="scnanoseq", metavar="string" , help="Project name for Seurat object."                                     ),
    make_option(c("-o", "--outdir"        ), type="character", default="./"       , metavar="path"   , help="Output directory."                                                   ),
    make_option(c("-r", "--outprefix"     ), type="character", default="seurat_qc", metavar="string" , help="Output prefix."                                                      )
)

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_matrix) && is.null(opt$input_dir)) {
    print_help(opt_parser)
    stop("Please provide either a single-cell/nuclei matrix or a directory containing a matrix.mtx, genes.tsv (or features.tsv) and barcodes.tsv.", call. = FALSE)
}

if (is.null(opt$flagstat)) {
    print_help(opt_parser)
    stop("Please provide the samtools flagstat file.", call. = FALSE)
}

######################
### READ IN INPUTs ###
######################

# Create the Seurat object
#NOTE: we do not perform any pre-filtering at this point

# cell or nuclei matrix (calling it cell for simplicity)

if (!is.null(opt$input_dir)) {
    cell_bc_matrix <- Read10X(data.dir = opt$input_dir,
                                gene.column = 1,
                                cell.column = 2)
    seurat_obj <- CreateSeuratObject(counts = cell_bc_matrix,
                                        min.cells = 0,
                                        min.features = 0,
                                        project = opt$id)
} else {
    cell_bc_matrix <- read.table(opt$input_matrix, sep="\t", header = TRUE, row.names = 1)
    seurat_obj <- CreateSeuratObject(counts = cell_bc_matrix,
                                        min.cells = 0,
                                        min.features = 0,
                                        project = opt$id)
}

flagstat_lines <- readLines(opt$flagstat)

#######################
### GET TOTAL READS ###
#######################

# This will grep the lines to get the line index containing the total read count
index_nums <- grep("total", flagstat_lines)

# This will parse out the total read count
total_reads <- as.numeric(gsub("([0-9]+).*$", "\\1", flagstat_lines[index_nums]))


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

# Generate the violin plot
violin_plot <- VlnPlot(object = seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), assay = "RNA")

# Generate the density plots
densityplot_nFeature <- plotSingleCellDensity(input_obj = seurat_obj, metadata_variable = "nFeature_RNA", scale_x_axis = TRUE)

densityplot_nCount <- plotSingleCellDensity(input_obj = seurat_obj, metadata_variable = "nCount_RNA", scale_x_axis = TRUE)

# Feature Scatter plot
feature_plot <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

# gather all plots
all_plots <- (violin_plot + densityplot_nFeature + densityplot_nCount + feature_plot)

# save plots
png(out_sample_qc_figs, width = 900, height = 900)

plot(all_plots)

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

# Get the "Median Feautures per Cell"
# Feature may be genes or transcripts
median_features_per_cell <- median(seurat_obj$nFeature_RNA)

# Get the "Total Number of Features" detected
total_number_features <- nrow(seurat_obj@assays$RNA@counts)

####################
### OUTPUT TABLE ###
####################

# Acquire all stats
output_table <- data.frame(est_cell_number, mean_reads_per_cell, median_features_per_cell, total_number_features)

# Set the colnames for the output table
colnames(output_table) <- c("Estimated Cell Number", "Mean Reads Per Cell", "Median Feautures Per Cell", "Total Number of Features")

# Write out the data frame
out_stats <- paste0(opt$outprefix,".csv")

write.csv(output_table, out_stats, row.names = FALSE)
saveRDS(seurat_obj, file = paste0(opt$outprefix,"_seurat.rds"))

####################
### Session Info ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
