# General Functions



################################################################################################
# Loads all required libraries
################################################################################################
load_libraries <- function () {
  # libs
  required_libs = c(
    c("Seurat", "readxl", "varhandle", "MASS", "dplyr", "tidyverse"),
    if ("decontaminate" %in% config$process) c("SoupX", "celda", "Matrix", "qlcMatrix", "FastCAR") else NULL,
    if ("analyse" %in% config$process) c("ggplot2", "xlsx", "reshape2","plotly","cowplot","patchwork") else NULL,
	if ("summarise" %in% config$process) c("ggplot2", "xlsx", "reshape2", "stringr", "ggforce") else NULL
  )
  required_libs = unique(required_libs)
		
  uninstalled = c()
  
  # looping through each lib, check if installed then loading
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
  }
}



################################################################################################
# Returns config list containing all required variables
  # Throws an error if a required config variable is empty
################################################################################################
get_config <- function(args, file=F) {
  if (file == F)
    config <- config::get(config = (if (length(args) > 0) args[[1]]), file = "config.yml")
  else
    config <- config::get(config = (if (length(args) > 0) args[[1]]), file=file)
  
  # checking for empty variables
  required = c("alpha", "threads", "quiet", "genes_ct_dotplots", "ct_order_dotplots", "pie_plot_cts", "input_dir",
               "output_dir", "process", "methods", "summary_histogram_labels", "sample_ids", "dataset")
                                    
  for (i in required) {
    if (config[[i]] == "" || is.null(config[[i]]))
      stop(paste("Config is empty @ ", i,sep=""))
  }

  
  # splitting "arrays"
  for (i in c("sample_ids", "genes_ct_dotplots", "ct_order_dotplots", "pie_plot_cts")) {
    config[[i]] = strsplit(config[[i]], split="|",fixed=T)[[1]]
  }
  
  # original cell annotations file is xlsx
  config$is_xlsx = sapply(config$methods, function(x) {
    return(config$recluster == F || x == "none" || x == "no_decontamination")
  })
  
  return(config)
}



################################################################################################
# Gets all file locations based on the config file and the current method
################################################################################################
get_files <- function (config, current_method, subset_index = F) {
  if (config$dataset == "mouse_kidney")
      repeat_names = config$sample_ids
    else if (config$dataset == "hgmm12k")
      repeat_names = c("hgmm12k")

  if (subset_index != F) {
    files = list(
      "CellRanger"        = paste(config$input_dir, subset_index, sep="/"),
      "CellRangerMerged"  = paste(config$input_dir, subset_index, sep="/"),
      "CellBender"        = paste(config$input_dir, "cellbender", subset_index, sep="/"),
      "output"            = paste(config$output_dir, current_method, subset_index, sep="/"),
      "CellAnnotations"   = paste(config$input_dir, config$original_cell_annotations, sep="/")
    )
  } else { 
    files = list(
      "CellRanger" = paste(config$input_dir, "CellRanger", sep="/"),
      "CellRangerMerged" = paste(config$input_dir, "CellRanger_merged/raw_gene_bc_matrices", sep="/"),
      "Filtered" = paste(config$input_dir, "Filtered_Feature_Barcode_Matrices", sep="/"),
      "CellBender" = paste(config$input_dir, "cellbender", sep="/"),
      "GeneSignatures" = paste(config$input_dir, config$gene_signatures, sep="/"),
      "output" = paste(config$output_dir, current_method, sep="/"),
      "OcraRelDir" = paste(config$ocra_dir, current_method, sep="/")
    )
    
    if (config$dataset == "mouse_kidney") {
      files[["CellAnnotations"]] = if (config$is_xlsx[[current_method]]) 
        paste(config$input_dir, config$original_cell_annotations, sep="/") 
      else 
        paste(config$input_dir, config$new_cell_annotations, sep="/")
    } else if (config$dataset == "hgmm12k") {
      files[["CellAnnotations"]] = paste(config$input_dir, config$original_cell_annotations, sep="/")
    }
    
    # adds filtered file paths for cellbender
    if (config$dataset == "mouse_kidney" & current_method == "cellbender") {
      files$Filtered = sapply(repeat_names, FUN = function(x) {
        paste(files$Filtered,"/", x,".txt", sep="")
      })
    }
    
    # specific locations based on the current method and the sample ids
    if (substring(current_method,0,5) == "soupx") {
      files$special = rep(files$GeneSignatures, length(config$sample_ids))
    } else if (current_method == "cellbender") {
      files$special = sapply(repeat_names, function(x) {
        paste(files$Filtered,"/", x,".txt", sep="")
      }, USE.NAMES=F)
    } else if (current_method == "no_decontamination") { 
      files$special = sapply(config$sample_ids, function(x) {
        paste(files$Filtered,"/", x,".txt", sep="")
      }, USE.NAMES=F)
    }
  }

  files$dir = sapply(repeat_names, FUN = function(x) {
    if (substring(current_method,0,5) == "soupx") 
      paste(files$CellRanger, x, sep="/")
    else if (substring(current_method,0,7) == "decontx")
      paste(files$Filtered,"/", x,".txt", sep="")
    else if (current_method == "cellbender")
      paste(files$CellBender, "/", x, "/", x, "_filtered.h5", sep="")
    else
      paste(files$CellRanger, x, sep="/")
  }, USE.NAMES=F)
  
  return(files)
}



									
									
									
									
################################################################################################
# Log function
# TODO improve
    # log to file?
    # everything should go through this, maybe a level as well? ie warning, error, etc.
################################################################################################
log_print <- function (msg) {
  # TODO output file as well?
  if (!config$quiet)
    print(msg)
}

									   

################################################################################################
# Checks if samples.combined is loaded, if not -> loads it from the file
################################################################################################
load_rda <- function (samples.combined, file, load_no_decont = FALSE) {
  if (load_no_decont) 
    samples.combined <- readRDS(paste(config$output_dir, "no_decontamination", file, sep="/"))
  
  else if (is.null(samples.combined)) {
    print("Reading Rda from file...")
    samples.combined <- readRDS(paste(files$output, file, sep="/"))
  }
	
  return(samples.combined)
}




