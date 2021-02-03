# command line args
args = commandArgs(trailingOnly=TRUE)

# loading functions from separate scripts
source("../improved_scripts/scripts/general_functions.R")
source("../improved_scripts/scripts/decontamination_functions.R")
source("../improved_scripts/scripts/analysis_functions.R")

if (!requireNamespace('config', quietly=T))
  stop("R Package 'config' not installed.")
# DO NOT LOAD CONFIG - CLASHES W/ SEURAT MERGE

# getting config and loading libs
config <- get_config(args)
load_libraries()


################################################################################################
# Function to decontaminate samples, save the feature barcode matrix, and return a seurat object
################################################################################################
decontaminate_samples <- function (config, files, current_method) {
  samples = as.list(config$sample_ids)
  names(samples) <- config$sample_ids
  
  ## Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
  for (i in seq(length(config$sample_ids))) {
    sample_id = config$sample_ids[i]
    
    print(paste("Starting",sample_id))
    samples[[sample_id]] = get_sample(sample_id, files$dir[i], current_method, files$CellAnnotations, files$special[i], 
									  !(config$recluster == F || (current_method == "none" || current_method == "no_decontamination")))

    # ensuring formatting of cell barcodes is the same (across all analyses)
    samples[[sample_id]]$seurat = fix_barcodes(samples[[sample_id]]$seurat)
  }
  
  save_matrices(samples, paste(files$output, "/matrices/",sep=""))
  
  ### Combining seurat objects
  print("Combining seurat objects")
  
  # loops through each object in samples <- adds metadata to `seurat.decont` 
  samples.seurat <- lapply(samples, function(x) {
    x$seurat@meta.data$orig.ident = x$sample_id
    x$seurat@meta.data$preservation = if (length(str_split(x$sample_id, "_",simplify=TRUE))>2) "MeOH" else "fresh"
    x$seurat@meta.data$method = "decont"
    return(x$seurat)
  })
  
  # merging
  samples.combined <- merge(samples.seurat[[1]], samples.seurat[2:6], add.cell.ids = config$sample_ids)

  # factoring metadata
  samples.combined@meta.data$orig.ident = factor(samples.combined@meta.data$orig.ident)
  samples.combined@meta.data$preservation = factor(samples.combined@meta.data$preservation)
  samples.combined@meta.data$method = factor(samples.combined@meta.data$method)
  

  ### Processing prior to integration
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



################################################################################################
# Integrating samples via Seurat::IntegrateData
################################################################################################
integrate_samples <- function (config, files, samples.combined) {
  print("Finding integration anchors")
  samples.anchors <- FindIntegrationAnchors(object.list = samples.combined)
 
  print("Integrating")
  samples.combined <- IntegrateData(anchorset = samples.anchors)
  
  # If reclustered <- saves new cell annotations
  if (config$recluster)
    write.table(as.matrix(Idents(samples.combined)), paste(files$output, "/new_clus.tsv", sep=""), sep="\t")
  
  DefaultAssay(samples.combined) <- "integrated"
  
  #### Dimension Reduction
  print("Dimension reduction")
  samples.combined <- ScaleData(samples.combined, verbose = FALSE)
  samples.combined <- RunPCA(samples.combined, npcs = 30, verbose = FALSE)
  samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims=1:30)
  
  # adding metadata <- mostly for plotting the dotplots
  samples.combined <- adding_metadata(samples.combined)
  
  # saves the object as an rda
  saveRDS(samples.combined,paste(files$output, "Rda/integrated_rd.Rda", sep="/"))
  return(samples.combined)
}



################################################################################################
# Analysing samples; mostly plotting
################################################################################################
analyse_samples <- function (config, files, samples.combined) {
  # Differentially expressed genes
  analyse_DEGs()
  # UMAP
  analyse_UMAPs()
  # Reclustered plots / tables / etc.
  if (config$recluster)
    analyse_recluster()
}



################################################################################################
# Summarising samples
################################################################################################
summarise_samples <- function (config, files, samples.combined) {
  # TODO
}



################################################################################################
#start
################################################################################################
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



# Looping through methods
for (current_method in config$methods) {
  print(paste(rep("#",30),collapse=""))
  print(paste("Starting",current_method))
  print(paste(rep("#",30),collapse=""))
  
  files = get_files(config, current_method)
  samples.combined=NULL
	
  # Decontamination
  if ("decontaminate" %in% config$process) {
    print("Decontaminating")
    samples.combined = decontaminate_samples(config, files, current_method)
    print("Decontamination completed")
  }
  
  # Integration
  if ("integrate" %in% config$process) {
	samples.combined <- load_rda(samples.combined, "Rda/decontaminated_samples.Rda")
    
    print("Integrating")
    samples.combined = integrate_samples(config, files, samples.combined)
    print("Integration completed")
  }
  
  # Analysis
  if ("analyse" %in% config$process) {
	samples.combined <- load_rda(samples.combined, "Rda/integrated_rd.Rda")
    
    print("Analysing")
    samples.combined = analyse_samples(config, files, samples.combined)
    print("Analysis completed")
  }
  
  # Summary
  if ("summarise" %in% config$process) {
	samples.combined <- load_rda(samples.combined, "Rda/integrated_rd.Rda")
    
    print("Summarising")
    samples.combined = summarise_samples(config, files, samples.combined)
    print("Summary completed")
  }
}
