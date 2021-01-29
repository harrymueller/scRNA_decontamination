#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading general functions
source("../improved_scripts/scripts/general_functions.R")
source("../improved_scripts/scripts/decontamination_functions.R")

load_libraries()

# getting config
config <- get_config(args)




decontaminate_samples <- function (config, files, current_method) {
  samples = as.list(config$sample_ids)
  names(samples) <- config$sample_ids
  
  ## Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
  #method = ["none", "soupx:autoEstCont", "soupx:gene"]
  for (i in seq(length(config$sample_ids))) {
    sample_id = config$sample_ids[i]
    
    print(paste("Starting",sample_id))
    samples[[sample_id]] = get_sample(sample_id, files$dir[i], current_method, files$CellAnnotations, files$special, config$recluster)
    
    # ensuring formatting of cell barcodes is the same (across all analyses)
    samples[[sample_id]]$seurat.decont = fix_barcodes(samples[[sample_id]]$seurat.decont, paste(output, "errors.txt", sep="/"))
  }
  
  save_matrices(samples, paste(files$output, "/matrices/",sep=""))
  save_Rda(samples, paste(files$output, "/Rda/decontaminated_samples/", ))
}

integrate_samples <- function (config) {
	
}

analyse_samples <- function (config) {

}

################
#start
################
# multithreading
if (config$threads != 1) {
  options(future.globals.maxSize = 4 * 1024^3)
  
  # Without \/, Seurat will produce errors based on randomisation of processes. Can be ignored (https://github.com/satijalab/seurat/issues/3622)
  options(future.rng.onMisuse="ignore")
  
  # Multiprocessing
  library(future)
  library(future.apply)
  plan("multiprocess", workers = config$threads)
}

for (current_method in config$methods) {
  files = get_files(config, current_method)
  if ("decontaminate" %in% config$process) {
    decontaminate(config, files, current_method)
  } 
  
}
# for each method...
# decont, integ, analyse ...
