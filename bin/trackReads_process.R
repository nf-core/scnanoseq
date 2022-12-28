#!/usr/bin/env Rscript

######################
### LOAD LIBRARIES ###
######################
library(optparse)
library(tidyverse)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
  make_option(c("-i", "--input_bambu_obj"), type="character", default=NULL            , metavar="path"   , help="Path to bambu R object containing readToTranscriptMaps metadata."),
  make_option(c("-g", "--gene_meta"      ), type="character", default=NULL            , metavar="path"   , help="Path to bambu's counts_transcript output."), # for gene ids
  make_option(c("-o", "--outdir"         ), type="character", default="./"            , metavar="path"   , help="Output directory."),
  make_option(c("-r", "--outprefix"      ), type="character", default="bambu_read_map", metavar="string" , help="Output prefix.")
)

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_bambu_obj)) {
  print_help(opt_parser)
  stop("Please provide a bambu R object.", call. = FALSE)
}

if (is.null(opt$gene_meta)) {
  print_help(opt_parser)
  stop("Please provide the bambu counts_transcript file.", call. = FALSE)
}

#######################
### CUSTOM FUNCTION ###
#######################

#TODO: bambu readtracks output does contain more than one row per read in some cases
# filter these upstream if NULL across and or not unique (leaving off for now for initial tests...)

#' process_trackreads
#' This function process the `readToTranscriptMaps` output from `bambu` to add annotations data
#' to each `readId`. This function is meant to annotate `equalMatches` and `compatibleMatches` separately.
#' The results can be merged to re-generate a single data.frame of all data.
#'
#' @param trackreads_input readToTranscriptMaps from bambu object - e.g.: `metadata(se)$readToTranscriptMaps`
#' @param annotation_index indices which match the annotations from bambu
#' @param variable_name variable name to process - either "equalMatches" or "compatibleMatches"
#'
#' @return A list of annotated `readToTranscriptMaps`
#' @export
#'
#' @examples
#' \dontrun{
#' process_trackreads(trackreads_input = metadata(se)$readToTranscriptMaps,
#'                    annotation_index = annotation_index,
#'                    variable_name = "equalMatches")
#' }
process_trackreads <- function(trackreads_input,
                               annotation_index,
                               variable_name = c("equalMatches", "compatibleMatches")) {
  # evaluate variable names
  variable_name <- match.arg(variable_name)
  
  # removes remove any spaces and `c()` to separate rows,
  # select the variable from input list and separate_rows
  
  trackReads <- lapply(X = trackreads_input, FUN = function(x) {
    
    trackReads <- as.data.frame(apply(x,2,as.character))
    
    # gsub any spaces, and `c()` from index column
    trackReads[[variable_name]] <- gsub(" ", "",
                                        gsub("[()]", "",
                                             gsub("c", "", trackReads[[variable_name]])
                                        )
    )
    
    # select, separate
    
    trackReads <- trackReads %>%
      dplyr::select(readId, {{ variable_name }}) %>%
      tidyr::separate_rows({{ variable_name }}, convert = TRUE)
    
    return(trackReads)
    
  })
  
  ### Add annotations based of index ###
  
  # make column names from unique variable names
  tx_name <- paste0(variable_name, "_TXNAME")
  gene_name <- paste0(variable_name, "_GENEID")
  
  index_IDs <- annotation_index %>%
    dplyr::rename({{ tx_name }} := TXNAME,
                  {{ gene_name }} := GENEID)
  
  index_IDs[[variable_name]] <- rownames(index_IDs)
  
  # merge to add annotations
  trackReads <- lapply(X = trackReads, FUN = function(x) {
    
    trackReads <- merge(x = x,
                        y = index_IDs,
                        by = variable_name,
                        all.x = TRUE,
                        sort = FALSE)
    
    return(trackReads)
    
  })
  
  return(trackReads)
  
}


######################
### READ IN INPUTs ###
######################

# cell or nuclei matrix (calling it cell for simplicity)
se <- readRDS(opt$input_bambu_obj)

annotation_index <- read.table(opt$gene_meta, header = TRUE)

##########################
### PREPARE ANNOTATION ###
##########################
# NOTE: indices match the annotations found in rowRanges(se)
# while they match outputs from txt file, for clarity we use
# rowRanges(se) as the basis of annotation below and bring gene ids from txt out

rowRanges_index <- as.data.frame(names(rowRanges(se)))

colnames(rowRanges_index) <- "TXNAME"

annotation_index <- dplyr::select(annotation_index, TXNAME, GENEID)

# sanity check that index order matches original rowRanges(se)
stopifnot(identical(rowRanges_index$TXNAME, annotation_index$TXNAME))

##########################
### PROCESS READTRACKS ###
##########################

# initially, process equalMatches and compatibleMatches separately
trackReads_list_equal <- process_trackreads(trackreads_input = metadata(se)$readToTranscriptMaps,
                                            annotation_index = annotation_index,
                                            variable_name = "equalMatches")

trackReads_list_compatible <- process_trackreads(trackreads_input = metadata(se)$readToTranscriptMaps,
                                                 annotation_index = annotation_index,
                                                 variable_name = "compatibleMatches")

###########################
### RE-GROUP BY READ ID ###
###########################
# NOTE: this may be a section (of re-grouping) that could need further optimization if needed

# one row per read name

# a little cleaner but takes longer
# obj %>% 
#   dplyr::group_by(readId) %>%
#   dplyr::summarise_all(list(paste = ~ paste(na.omit(.), collapse = ",")))

trackReads_list_equal <- lapply(X = trackReads_list_equal, FUN = function(x) {
  
  trackReads_grouped <- x %>%
    dplyr::group_by(readId) %>%
    dplyr::mutate(equalMatches = paste(equalMatches, collapse = ","),
                  equalMatches_TXNAME = paste(equalMatches_TXNAME, collapse = ","),
                  equalMatches_GENEID = paste(equalMatches_GENEID, collapse = ",")
    )
  
  trackReads_grouped <- unique(trackReads_grouped)
  
  return(trackReads_grouped)
  
})


trackReads_list_compatible <- lapply(X = trackReads_list_compatible, FUN = function(x) {
  
  trackReads_grouped <- x %>%
    dplyr::group_by(readId) %>%
    dplyr::mutate(compatibleMatches = paste(compatibleMatches, collapse = ","),
                  compatibleMatches_TXNAME = paste(compatibleMatches_TXNAME, collapse = ","),
                  compatibleMatches_GENEID = paste(compatibleMatches_GENEID, collapse = ",")
    )
  
  trackReads_grouped <- unique(trackReads_grouped)
  
  return(trackReads_grouped)
  
})

########################################################
### MERGE EQUAL AND COMPATIBLE MATCHES BACK TOGETHER ###
########################################################

trackReads_list_all <- mapply(FUN = function(x, y) {
  
  trackReads <- merge(x = x,
                      y = y,
                      by = "readId",
                      all = TRUE,
                      sort = FALSE)
  
  return(trackReads)
  
}, x=trackReads_list_equal, y=trackReads_list_compatible, SIMPLIFY = FALSE)

###################################
### MATCH READ ORDER FROM INPUT ###
###################################

# not a requirement, but for consistency re-ordering based off original read id order
# NOTE: some reads in bambu output are not in the same order (e.g some appear later on
# than prior occurrences). Thus this outputs are ordered by order of appearance of readId

trackReads_list_all <- mapply(FUN = function(x, y) {
  
  unique_read_order <- unique(y$readId)
  
  # safety check that all reads are present in both sets of data:
  stopifnot(all(unique_read_order %in% x$readId))
  
  # re-order
  trackReads <- x[match(unique_read_order, x$readId),]
  
  return(trackReads)
  
}, x=trackReads_list_all, y=metadata(se)$readToTranscriptMaps, SIMPLIFY = FALSE)


####################
### SAVE OUTPUTS ###
####################

# setting work dir to output dir
if (file.exists(opt$outdir) == FALSE) {
  dir.create(opt$outdir, recursive=TRUE)
}

setwd(opt$outdir)

# save annotated object with all data
saveRDS(trackReads_list_all, file = "./scnanoseq_trackReads_process.rds")

# save each element of the list for per sample outputs
# save all and transcript / gene id only
mapply(FUN = function(x,y) {
  
  file_name_all <- paste0("./",y,"_trackreads_all.txt")
  file_name_txname <- paste0("./",y,"_trackreads_txname.txt")
  file_name_geneid <- paste0("./",y,"_trackreads_geneid.txt")
  
  tx_data <- dplyr::select(x, readId,
                           equalMatches, equalMatches_TXNAME,
                           compatibleMatches, compatibleMatches_TXNAME)
  
  gene_data <- dplyr::select(x, readId,
                             equalMatches, equalMatches_GENEID,
                             compatibleMatches, compatibleMatches_GENEID)
  
  # save all
  write.table(x,
              file_name_all, sep = "\t",
              row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # save transcripts only
  write.table(tx_data,
              file_name_txname, sep = "\t",
              row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # save genes only
  write.table(gene_data,
              file_name_geneid, sep = "\t",
              row.names = FALSE, quote = FALSE, col.names = TRUE)
  
}, x = trackReads_list_all, y = names(trackReads_list_all))


####################
### Session Info ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
