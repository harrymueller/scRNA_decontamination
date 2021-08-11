# source files
source("scripts/general_functions.R")
source("scripts/mouse_kidney/analysis_functions.R")

# template config and files
config <- get_config("subset_reclus_pt_cells")

# swap these two due to the way `load_rda` works
output = config$output_dir
config$output_dir = config$input_dir

load_libraries()

degs_small_clus <- function(pt_subset, method, p, id) {
  DefaultAssay(pt_subset) = "RNA"
  # canabilised from scripts/mouse_kidney/DEGs.R
  pt_subset@meta.data["small_clus"] = factor(pt_subset$seurat_clusters == id, levels = c(T, F), labels=c("small_clus", "not_small_clus"))
  Idents(pt_subset) <- "small_clus"
  print(xtabs(~Idents(pt_subset)))
  return
  m = FindMarkers(pt_subset, ident.1 = "small_clus", ident.2 = "not_small_clus", 
                  verbose = F, logfc.threshold = 1, min.pct=0.5,test.use = "wilcox")
  m = m[m$p_val_adj < 0.05,]
  m = data.frame("names" = rownames(m), m[,1:5])
  m = m[order(m$avg_logFC),]
  # save to file
  wb <- createWorkbook()
  
  s <- createSheet(wb, "overexpressed_genes_>8_celltypes")
  addDataFrame(c(
    paste("DEGs between cluster", id, "for", p, "cells")
  ), s, row.names = FALSE, col.names = FALSE)
  
  # DEGs by cell type
  addDataFrame(m, s)
  
  name = paste("degs_", method, "_",p, ".xlsx", sep="")
  saveWorkbook(wb, paste(output, method, name, sep="/"))
  DefaultAssay(pt_subset) = "integrated"
}

dotplot <- function(sample) {
  DefaultAssay(sample) = "RNA"
  features = c("Kap", "Glyat", "Gpx3", "Gpx1", "Emcn", "Flt1", "Ly6c1", "Ifitm3")
  s=20
  plot = DotPlot(sample, features = features, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("CellBender: Expression of Genes in Seurat Clusters of PT Cells") + theme(text=element_text(size=12, family="TT Times New Roman"))

  name = paste("gene expr dotplot.png",sep="")
  ggsave(filename = paste(output, m, name, sep="/"), plot, width=9, height=7)  
  DefaultAssay(sample) = "integrated"
}

# loop through decont methods inc. no decont
for (m in config$methods) {
  files = get_files(config, m)
  
  # load data and add metadata
  samples.combined <- load_rda(NULL, "Rda/integrated_rd.Rda")
  samples.combined <- adding_metadata(samples.combined)

  Idents(samples.combined) = samples.combined$celltype
  DefaultAssay(samples.combined) = "integrated"
  
  # make dir for plots
  if (!dir.exists(paste(output, m, sep="/")))
    dir.create(paste(output, m, sep="/"))
  
  # seperately for fresh and meoh fixed
  for (p in c("fresh", "MeOH", "all")) {
  #for (p in c("all")) {
    pt_subset = samples.combined[,samples.combined$celltype == "PT"]
    
    # seperate by pres type then select only pt cells
    if (p != "all")
      pt_subset = pt_subset[,pt_subset$preservation == p]

    # recluster
    pt_subset <- FindNeighbors(pt_subset, dims = 1:10)
    pt_subset <- FindClusters(pt_subset, resolution = 0.5)
    
    # DEGs for small cluster for cellbender
    if (m == "cellbender" && p != "fresh") 
      degs_small_clus(pt_subset, m, p, 10)
    
    if (m == "cellbender" && p == "all")
      dotplot(pt_subset)
    
    Idents(pt_subset) = pt_subset$seurat_clusters
    
    pt_subset <- RunUMAP(pt_subset, reduction = "pca", dims=1:30)
    
    # umap & save
    plot = DimPlot(pt_subset, reduction = "umap",label=F) + 
      theme(text=element_text(size=16, family="TT Times New Roman")) +
      ggtitle(paste(m, ": After subsetting and reclustering (", p, ")", sep=""))
    
    name = paste(m, "_", p, "_plot.png",sep="")
    ggsave(filename = paste(output, m, name, sep="/"), plot, width=9, height=7)
  }
}



