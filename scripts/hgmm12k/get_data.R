################################################################################################
# Returns a seurat object containing the raw hgmm12k data
################################################################################################
get_raw_hgmm <- function(dir, types) {
  hg = Read10X(paste(dir, "raw_gene_bc_matrices", types[1], sep="/"))
  mm = Read10X(paste(dir, "raw_gene_bc_matrices", types[2], sep="/"))
  
  combined = rbind(hg, mm)
  return(CreateSeuratObject(combined))
}

################################################################################################
# Returns a seurat object containing the filtered hgmm12k data
################################################################################################
get_filtered_hgmm <- function(dir, annotation_path, types) {
  # raw data
  hg = Read10X(paste(dir, "raw_gene_bc_matrices", types[1], sep="/"))
  mm = Read10X(paste(dir, "raw_gene_bc_matrices", types[2], sep="/"))
  
  gem_classifications = read.csv(annotation_path)
  
  # using gem classifications to filter raw data - also removes multiplets
  hg.f = hg[,colnames(hg) %in% gem_classifications$barcode[gem_classifications$call != "Multiplet"]]
  mm.f = mm[,colnames(mm) %in% gem_classifications$barcode[gem_classifications$call != "Multiplet"]]
  
  combined = rbind(hg.f, mm.f)

  # removing genes w/ less than 10 counts across all filtered cells
  combined = combined[names(which(Matrix::rowSums(combined) > 10)),]
  
  # setting idents
  cell_annotations <- c(unfactor(gem_classifications$call[gem_classifications$call != "Multiplet"]))
  combined <- CreateSeuratObject(combined)
  Idents(combined) <- cell_annotations
  
  return(combined)
}
