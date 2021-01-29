# GENERAL FUNCTIONS
load_libraries <- function () {
  required_libs = c(
    c("Seurat", "readxl", "varhandle", "MASS", "dplyr", "tidyverse"),
    if ("decontaminate" %in% config$process) c("SoupX", "celda") else NULL,
    if ("analyse" %in% config$process) c("ggplot2", "xlsx", "reshape2","plotly","cowplot","patchwork") else NULL
  )
  uninstalled = c()
  
  for (i in required_libs) {
    if (!requireNamespace(i, quietly=T)) {
      # lib not found
      uninstalled = append(uninstalled, i)
    } else {
      # loading lib
      suppressMessages({suppressWarnings({
        library(i, character.only=T)
      })})
    }
  }
  
  if (!is.null(uninstalled)) {
    stop(paste("The following packages are required but not installed:",paste(uninstalled, collapse=", ")))
  } else {
    config = data.frame(quiet=F)
    log_print("All libraries loaded successfully")
  }
}

get_config <- function(args, testing=F) {
  if (!testing)
    config <- config::get(config = (if (length(args) > 0) args[[1]]))
  else
    config <- config::get(config="harry")
  
  # TODO check for empty strings in critical places
  
  # TODO split rest of "arrays"
  config$sample_ids = strsplit(config$sample_ids, split = "|",fixed=T)[[1]]
  
  return(config)
}

log_print <- function (msg) {
  # TODO output file as well?
  if (!config$quiet)
    print(msg)
}



get_files <- function (config, current_method) {
  files = list(
    "CellRanger" = paste(config$input_dir, "CellRanger", sep="/"),
    "Filtered" = paste(config$input_dir, "Filtered_Barcode_Matrices", sep="/"),
    "CellBender" = paste(config$input_dir, "CellBender", sep="/"),
    "CellAnnotations" = if (current_method == "none") paste(config$input_dir, config$original_cell_annotations_file, sep="/") 
                        else paste(config$input_dir, config$new_cell_annotations_file, sep="/"),
    "GeneSignatures" = paste(config$input_dir, "gene_signatures.xlsx", sep="/"),
    "output" = paste(config$output_dir, current_method, sep="/")
  )
  
  if (substring(current_method,0,5) == "soupx") {
    files$dir = sapply(config$sample_ids, FUN = function(x) {
      paste(files$CellRanger, x, sep="/")
    }, USE.NAMES=F)
  } else if (substring(current_method,0,7) == "decontx") {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  } else if (current_method == "cellbender") {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$CellBender, "/", x, "/", x, "_filtered.h5", sep="")
    }, USE.NAMES=F)
    
    files$special = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  } else {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$CellRanger, x, sep="/")
    }, USE.NAMES=F)
    
    files$special = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  }
  
  return(files)
}

files <- get_files(config, "soupx")


