# source files
source("scripts/general_functions.R")

# template config and files
config <- get_config("subset_reclus_pt_cells")

# swap these two due to the way `load_rda` works
output = config$output_dir
config$output_dir = config$input_dir

load_libraries()

# loop through decont methods inc. no decont
for (m in config$methods) {
  files = get_files(config, m)
  
  # load data and add metadata
  samples.combined <- load_rda(NULL, "Rda/integrated_rd.Rda")
  samples.combined <- adding_metadata(samples.combined)
  
  # ignoring meoh vs fresh
  Idents(samples.combined) = samples.combined$celltype
  DefaultAssay(samples.combined) = "integrated"
  # select only PT cells
  pt_subset = samples.combined[,samples.combined$celltype == "PT"]
  
  # recluster
  pt_subset <- FindNeighbors(pt_subset, dims = 1:10)
  pt_subset <- FindClusters(pt_subset, resolution = 0.5)
  Idents(pt_subset) = pt_subset$seurat_clusters
  
  pt_subset <- RunUMAP(pt_subset, reduction = "pca", dims=1:30)
  
  # umap & save
  p = DimPlot(pt_subset, reduction = "umap",label=F) + 
    theme(text=element_text(size=16, family="TT Times New Roman")) +
    ggtitle(paste(m, ": After subsetting and reclustering", sep=""))
  
  ggsave(paste(output, "/", m, "_plot.png",sep=""),p,width=9,height=7)
}


