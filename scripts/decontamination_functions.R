get_sample <- function(sample_id, dir, method, cell_annotations_path, special_files, use_new_clus) {
  ## Cell Annotations
  if (method != "none")
    cell_annotations = get_clusters(cell_annotations_path, sample_id, use_new_clus)
  
  out = list("sample_id"=sample_id)
  
  # SOUPX
  if (substring(method,0,5)=="soupx") {
    # loads dir in 'SoupChannel' object
    sc = load10X(dir)
    sc = setClusters(sc, setNames(cell_annotations, names(cell_annotations)))
    
    # removing any genes in the SC that are not in the filtered data
    if (length(sc$soupProfile[,1]) != length(rownames(sc$toc))) {
      sc$soupProfile = sc$soupProfile[rownames(sc$toc),]
    }
    
    if (method == "soupx:autoEstCont") {
      sc = autoEstCont(sc)
    } else if (method =="soupx:background_genes") {
      markers = get_markers(special_files, 3)
      useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = markers)
      sc = calculateContaminationFraction(sc,markers,useToEst=useToEst, forceAccept = T)
    } else if (method == "soupx:top_background_genes") {
      if (special_files == FALSE) {
        stop("No special_files given, but path to gene_signatures is required for soupX")
      }
      
      markers = get_top_n_markers(special_files, sc, 25)
      useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = markers)
      sc = calculateContaminationFraction(sc,markers,useToEst=useToEst, forceAccept = T)
    }
    
    decont_matrix = adjustCounts(sc)
    cont_matrix = sc$toc
  } else if (substring(method,0,7) == "decontx") {
    filtered = read.csv(dir,header = TRUE,sep = "\t")
    cont_matrix = as.matrix(filtered)
    
    if (method == "decontx:with_cell_types") {
      ## Cell Annotations
      cell_annotations = as.numeric(factor(cell_annotations))
      
      decont_matrix = decontX(cont_matrix, z=cell_annotations)$resList$estNativeCounts
    } else {
      decont_matrix = decontX(cont_matrix)$resList$estNativeCounts
    }
    
  } else if (method == "cellbender") {
    if (special_files == FALSE) {
      stop("No special_files given, but path to filtered files is required for no decontamination")
    }
    
    filtered = read.csv(special_files,header = TRUE,sep = "\t")
    cont_matrix = as.matrix(filtered)
    
    decont_matrix <- Read10X_h5(dir,use.names=T)
    
    # formatting cell barcodes to be '<sample id>_<barcode>'
    dimnames(decont_matrix)[[2]] = sapply(dimnames(decont_matrix)[[2]], function(x) {paste(sample_id, substring(x, 0, nchar(x)-2), sep="_")}, USE.NAMES=F)
    
    # only keep filtered cells and genes
    decont_matrix = decont_matrix[,(dimnames(decont_matrix)[[2]] %in% dimnames(cont_matrix)[[2]])]
    decont_matrix = decont_matrix[(dimnames(decont_matrix)[[1]] %in% dimnames(cont_matrix)[[1]]),]
    
    # reformatting cell barcodes to match barcodes in cell_annotations
    l = sapply(str_split(colnames(decont_matrix),"_"), function(n) paste(tail(n,1),"-1",sep=""))
    colnames(decont_matrix) = l
    
    # "fixing" cell annotations
    cell_annotations = cell_annotations[names(cell_annotations) %in% names(Idents(out$seurat.decont))]
    cell_annotations = cell_annotations[order(match(names(cell_annotations), names(Idents(out$seurat.decont))))]
  } else {
    filtered = read.csv(special_files,header = TRUE,sep = "\t")
    decont_matrix = as.matrix(filtered) # not actually "decontaminated" - however named this way
  } 
  
  out$seurat.decont = CreateSeuratObject(decont_matrix)
  
  # Add cluster information
  if (method != "none")
    Idents(out$seurat.decont) <- cell_annotations
  
  return(out)
}