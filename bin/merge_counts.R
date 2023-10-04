#!/usr/bin/env Rscript

################
### PACKAGES ###
################

library(optparse)
library(tidyverse)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(c("-i", "--input_mtx_1"), type="character", default=NULL, metavar="path", help="Input matrix to merge"),
    make_option(c("-j", "--input_mtx_2"), type="character", default=NULL, metavar="path", help="Input matrix to merge"),
    make_option(c("-o", "--output"      ), type="character", default=NULL, metavar="path", help="Output matrix"))

opt_parser <- OptionParser(option_list = params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_mtx_1)){
    print_help(opt_parser)
    stop("Please provide input matrix 1")
}

if (is.null(opt$input_mtx_2)){
    print_help(opt_parser)
    stop("Please provide input matrix 2")
}

if (is.null(opt$output)){
    print_help(opt_parser)
    stop("Please provide an output file")
}


####################
### DO THE THING ###
####################

# Read the files into tables
message("--- Convert inputs to tables ---")

input_table_1 <- read.table(opt$input_mtx_1, sep="\t", header = TRUE)
input_table_2 <- read.table(opt$input_mtx_2, sep="\t", header = TRUE)

# Full join the table on the gene column
message("--- Full join the table ---")

merged_table <- full_join(input_table_1, input_table_2, by = "gene", suffix = c(".x", ".y"))

# Make sure to replace the NAs with 0s
merged_table <- merged_table %>% replace(is.na(.), 0)

# Get the list of 'correct' columns, i.e. don't have the '.x' or '.y' suffix
message("--- Get the correct list of columns ---")

# For the summarization, remove the column containing strings
message("--- Remove the columns ---")

shared_column_table <- merged_table %>% select(ends_with(".x") | ends_with(".y"))
unshared_column_table <- merged_table %>% select(!(ends_with(".x") | ends_with(".y")))

# We only want to look at merged
bc_prefixes <- unique(sub("\\..*", "", colnames(shared_column_table)))

# Do rowsums for columns that are the same barcode
summed_mtx <- sapply(bc_prefixes, function(x) rowSums(shared_column_table[,startsWith(colnames(shared_column_table), x)]))

# Lets merge the columns back together
merged_mtx <- cbind(unshared_column_table, summed_mtx)

##############
### OUTPUT ###
##############

# Write out the data frame
message("--- Output the table ---")
write.table(merged_mtx, opt$output, sep = '\t', row.names = FALSE, quote = FALSE)

####################
### Session Info ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()

