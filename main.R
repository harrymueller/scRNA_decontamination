# command line args
args = commandArgs(trailingOnly=TRUE)

# loading functions from separate scripts
source("scripts/general_functions.R")
source("scripts/decontamination_functions.R")
source("scripts/clustering.R")

# checking for 'config' pacakge
if (!requireNamespace('config', quietly=T))
  stop("R Package 'config' not installed.") # DO NOT LOAD CONFIG - CLASHES W/ SEURAT MERGE

# getting config and loading libs
print("Getting the config and loading libraries...")
config <- get_config(args)
load_libraries()
print("Config and libraries loaded successfully.")

# loading dataset-specific functions
if (config$dataset == "mouse_kidney") {
  source("scripts/mouse_kidney/analysis_functions.R")
  source("scripts/mouse_kidney/DEGs.R")
  source("scripts/mouse_kidney/summarising_functions.R")
} else if (config$dataset == "hgmm12k") {
  source("scripts/hgmm12k/get_data.R")
}

################################################################################################
# Function to decontaminate samples, save the feature barcode matrix, and return a seurat object
################################################################################################
decontaminate_samples <- function (current_method) {
  # checks if folders are created - if not makes them
  paths = sapply(c("", "Rda", "Rda/decontaminated_samples", "matrices"), 
                 function(x) paste(files$output,x, sep="/"), USE.NAMES=FALSE)
  for (p in paths) 
    if (!dir.exists(p))
      dir.create(p)
  
  # creating list to store decontaminated seurat objects
  if (config$dataset == "mouse_kidney") {
    samples = as.list(config$sample_ids)
    names(samples) <- config$sample_ids
    s = seq(length(config$sample_ids))
  } 
  else if (config$dataset == "hgmm12k") {
    samples = as.list("hgmm12k")
    names(samples) = "hgmm12k"
    s = seq(1)
  }
  
  ## Creates and saves individual R list objects 
  for (i in s) {
    sample_id = names(samples)[i]
    
    print(paste("Starting",sample_id))
    samples[[sample_id]] = get_sample(i, sample_id, current_method)
    
    # ensuring formatting of cell barcodes is the same (across all analyses)
    samples[[sample_id]]$seurat = fix_barcodes(samples[[sample_id]]$seurat)
  }
  
  save_matrices(samples)
  
  ### Combining seurat objects
  print("Combining seurat objects")
  
  # loops through each object in samples <- adds metadata to `seurat.decont` 
  samples.seurat <- lapply(samples, function(x) {
    x$seurat@meta.data$orig.ident = x$sample_id
    
    if (config$dataset == "mouse_kidney")
      x$seurat@meta.data$preservation = if (length(str_split(x$sample_id, "_",simplify=TRUE))>2) "MeOH" else "fresh"
        
    x$seurat@meta.data$method = "decont"
    return(x$seurat)
  })
  
  # merging seurat objects
  if (config$dataset == "mouse_kidney")
	  samples.combined <- merge(samples.seurat[[1]], samples.seurat[2:6], add.cell.ids = config$sample_ids)
  else
	  samples.combined <- samples.seurat[[1]]
	
  # factoring metadata
  samples.combined@meta.data$orig.ident = factor(samples.combined@meta.data$orig.ident)
  if (config$dataset == "mouse_kidney")
    samples.combined@meta.data$preservation = factor(samples.combined@meta.data$preservation)
  samples.combined@meta.data$method = factor(samples.combined@meta.data$method)
  

  ### Processing prior to integration
  print("Processing prior to integration")
  # splits combined seurat object into each individual seurat
  samples.combined = SplitObject(samples.combined, split.by="orig.ident")
  
  # normalises and finds variable features of individual seurats
  samples.combined <- lapply(X = samples.combined, FUN = function(x) {
    x <- NormalizeData(x)
    if (config$recluster) x <- reCluster(x, files$GeneSignatures, config$alpha)
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
integrate_samples <- function (samples.combined) {
  # Integrate samples
  if (config$dataset == "mouse_kidney") { # TODO <- make it dependant on length of samples.combined
    print("Finding integration anchors")
    samples.anchors <- FindIntegrationAnchors(object.list = samples.combined)

    print("Integrating")
    samples.combined <- IntegrateData(anchorset = samples.anchors)

    DefaultAssay(samples.combined) <- "integrated"
  } else
    samples.combined <- samples.combined$hgmm12k
  # If reclustered <- saves new cell annotations
  if (config$recluster)
    write.table(as.matrix(Idents(samples.combined)), paste(files$output, "/new_clus.tsv", sep=""), sep="\t")
  
  #### Dimension Reduction
  print("Dimension reduction")
  samples.combined <- ScaleData(samples.combined, verbose = FALSE)
  samples.combined <- RunPCA(samples.combined, npcs = 30, verbose = FALSE)
  samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims=1:30)
  
  # saves the object as an rda
  saveRDS(samples.combined,paste(files$output, "Rda/integrated_rd.Rda", sep="/"))
  return(samples.combined)
}



################################################################################################
# Analysing samples; mostly plotting
################################################################################################
analyse_samples <- function (samples.combined) {
  # checks for dir
  if (!dir.exists(paste(files$output, "plots", sep="/")))
    dir.create(paste(files$output, "plots", sep="/"))

  if (config$dataset == "mouse_kidney") # adding metadata
    samples.combined <- adding_metadata(samples.combined)
  # TODO FIX FOR mouse_kidney
  # UMAP
  Idents(samples.combined) = "celltype"
  p = DimPlot(samples.combined, reduction = "umap",label=F)
  ggsave(paste(files$output, "/plots/umap_plot.png",sep=""),p,width=9,height=7)

  # Mouse_Kidney analysis
  if (config$dataset == "mouse_kidney") {
    saveRDS(samples.combined, "~/Downloads/temp.Rda")
    # Differentially expressed genes
    analyse_DEGs(samples.combined)

    # Plotting pie charts of cell types
    plot_pie_ct(samples.combined, current_method, "preservation")
    plot_pie_ct(samples.combined, current_method, "method")

    # Reclustered plots / tables / etc.
    if (config$recluster == T)
      analyse_recluster(samples.combined, current_method)
  }
  
  # hgmm12k analysis
  else if (config$dataset == "hgmm12k") {
    # TODO ...
  }
}



################################################################################################
# Summarising samples
################################################################################################
summarise_samples <- function () {
  # summarise mouse_kidney analysis
  if (config$dataset == "mouse_kidney") {
    # check dir exists
    if (!dir.exists(paste(config$output_dir, "summary", sep="/")))
      dir.create(paste(config$output_dir, "summary", sep="/"))

    # summary histogram + summary degs
    deg_summary()

    if (config$recluster) {
      # ARI / NMI -> 1 doc & histograms
      ari_nmi = concat_ari_nmi()
      plot_ari_nmi()
    }
  }
  
  # summarise hgmm12k dataset
  else if (config$dataset == "hgmm12k") {
    # TODO ... 
  }
}



################################################################################################
# START
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


# If any methods in config$process...
if (length(intersect(c("decontaminate", "integrate", "analyse"), config$process)) > 0) {
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
      samples.combined = decontaminate_samples(current_method)
      print("Decontamination completed")
    }

    # Integration
    if ("integrate" %in% config$process) {
    samples.combined <- load_rda(samples.combined, "Rda/decontaminated_samples.Rda")

      print("Integrating")
      samples.combined = integrate_samples(samples.combined)
      print("Integration completed")
    }

    # Analysis
    if ("analyse" %in% config$process) {
      samples.combined <- load_rda(samples.combined, "Rda/integrated_rd.Rda")

      print("Analysing")
      samples.combined = analyse_samples(samples.combined)
      print("Analysis completed")
    }
  }
}
                 
# Summary
if ("summarise" %in% config$process) {    
  print("Summarising")
  samples.combined = summarise_samples()
  print("Summary completed")
}
