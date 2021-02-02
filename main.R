args = commandArgs(trailingOnly=TRUE)

# loading general functions
source("../improved_scripts/scripts/general_functions.R")
source("../improved_scripts/scripts/decontamination_functions.R")
source("../improved_scripts/scripts/analysis_functions.R")

if (!requireNamespace('config', quietly=T))
  stop("R Package 'config' not installed.")
# DO NOT LOAD CONFIG - CLASHES W/ SEURAT MERGE
config <- get_config(args)
if (F) {
#load_libraries()
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
  log_print("All libraries loaded successfully")
}
}
load_libraries()

decontaminate_samples <- function (config, files, current_method) {
  samples = as.list(config$sample_ids)
  names(samples) <- config$sample_ids
  
  ## Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
  #method = ["none", "soupx:autoEstCont", "soupx:gene"]
  for (i in seq(length(config$sample_ids))) {
    sample_id = config$sample_ids[i]
    
    print(paste("Starting",sample_id))
    samples[[sample_id]] = get_sample(sample_id, files$dir[i], current_method, files$CellAnnotations, files$special[i], !(config$recluster == F || (current_method == "none" || current_method == "no_decontamination")))
                                      # ^ Dont use new clusters if current method is no decontamination
    # ensuring formatting of cell barcodes is the same (across all analyses)
    samples[[sample_id]]$seurat = fix_barcodes(samples[[sample_id]]$seurat)
  }
  
  #save_matrices(samples, paste(files$output, "/matrices/",sep=""))
  
  ################################################################################################
  ### Combining seurat objects
  ################################################################################################
  print("Combining seurat objects")
  
  # loops through each object in samples <- adds metadata to `seurat.decont` 
  samples.seurat <- lapply(samples, function(x) {
    x$seurat@meta.data$orig.ident = x$sample_id
    x$seurat@meta.data$preservation = if (length(str_split(x$sample_id, "_",simplify=TRUE))>2) "MeOH" else "fresh"
    x$seurat@meta.data$method = "decont"
    return(x$seurat)
  })
  

  # combining
  samples.combined <- merge(samples.seurat[[1]], samples.seurat[2:6], add.cell.ids = config$sample_ids)

  # factoring metadata
  samples.combined@meta.data$orig.ident = factor(samples.combined@meta.data$orig.ident)
  samples.combined@meta.data$preservation = factor(samples.combined@meta.data$preservation)
  samples.combined@meta.data$method = factor(samples.combined@meta.data$method)
  
  ################################################################################################
  ### Processing prior to integration
  ################################################################################################
  print("Processing prior to integration")
  # splits combined seurat object into each individual seurat
  samples.combined = SplitObject(samples.combined, split.by="orig.ident")
  
  # normalises and finds variable features of individual seurats
  samples.combined <- lapply(X = samples.combined, FUN = function(x) {
    x <- NormalizeData(x)
    if (config$recluster) {
      x <- reCluster(x, files$GeneSignatures, config$alpha)
    }
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
  })
  
  saveRDS(samples.combined, paste(files$output, "Rda/decontaminated_samples.Rda", sep="/"))
  print("Completed decontamination")
  return(samples.combined)
}



# Integrating samples
integrate_samples <- function (config, files, samples.combined) {
  print("Finding integration anchors")
  samples.anchors <- FindIntegrationAnchors(object.list = samples.combined)
 
  print("Integrating")
  samples.combined <- IntegrateData(anchorset = samples.anchors)
  
  # Saving new cts
  if (config$recluster)
    write.table(as.matrix(Idents(samples.combined)), paste(files$output, "/new_clus.tsv", sep=""), sep="\t")
  
  DefaultAssay(samples.combined) <- "integrated"
  
  #### Dimension Reduction
  print("Dimension reduction")
  
  samples.combined <- ScaleData(samples.combined, verbose = FALSE)
  samples.combined <- RunPCA(samples.combined, npcs = 30, verbose = FALSE)
  samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims=1:30)
  
  samples.combined <- adding_metadata(samples.combined)
  
  saveRDS(samples.combined,paste(files$output, "Rda/integrated_rd.Rda", sep="/"))

  return(samples.combined)
}

analyse_samples <- function (config) {
  analyse_DEGs()
  analyse_UMAPs()
  if (config$recluster)
    analyse_recluster()
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
  print(paste(rep("#",30),collapse=""))
  print(paste("Starting",current_method))
  print(paste(rep("#",30),collapse=""))
  
  files = get_files(config, current_method)
  samples.combined=NULL
  if ("decontaminate" %in% config$process) {
    print("Decontaminating")
    samples.combined = decontaminate_samples(config, files, current_method)
    print("Decontamination completed")
  }
  
  if ("integrate" %in% config$process) {
    if (is.null(samples.combined)) {
      print("Reading Rda from file...")
      samples.combined <- readRDS(paste(files$output, "Rda/decontaminated_samples.Rda", sep="/"))
    }
    
    print("Integrating")
    samples.combined = integrate_samples(config, files, samples.combined)
    print("Integration completed")
  }
  
  if ("analyse" %in% config$process) {
    if (is.null(samples.combined)) {
      print("Reading Rda from file...")
      samples.combined <- readRDS(paste(files$output, "Rda/integrated_rd.Rda", sep="/"))
    }
    
    print("Analysing")
    samples.combined = analyse_samples(config, files, samples.combined)
    print("Analysis completed")
  }
  
}
# for each method...
# decont, integ, analyse ...
