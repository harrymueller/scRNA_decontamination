# All functions related to decontaminating samples



################################################################################################
# Returns a list containing a $seurat & $sample_id
################################################################################################
get_sample <- function(i, sample_id, method) {
  dir = files$dir[i]
  
  ## Cell Annotations
  if (method != "none")
    cell_annotations = get_clusters(files$CellAnnotations, sample_id, config$is_xlsx[[current_method]])
  
  out = list("sample_id"=sample_id)
  
  # SOUPX
  if (substring(method,0,5)=="soupx") {
    # loads dir in 'SoupChannel' object
    if (sample_id != "hgmm12k") {
      sc = load10X(dir)
    
      marker_file = files$special[i] # setting marker file as dir to gene_sig file
    }
    else if (sample_id == "hgmm12k") {
      filtered = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)
      marker_file = filtered # setting marker file as seurat object
      
      # converting seurat objects to sparse matrices to create SC obj
      filtered = filtered@assays$RNA@counts
      raw = get_raw_hgmm(files$CellRanger, config$sample_ids)@assays$RNA@counts
      sc = SoupChannel(raw, filtered)
    }
	  
    sc = setClusters(sc, setNames(cell_annotations, names(cell_annotations))) 
    
    # removing any genes in the SC that are not in the filtered data
    if (length(sc$soupProfile[,1]) != length(rownames(sc$toc))) {
      sc$soupProfile = sc$soupProfile[rownames(sc$toc),]
    }
    
    if (method == "soupx:autoEstCont") {
      sc = autoEstCont(sc)
    } else if (method =="soupx:background_genes" | method == "soupx:top_background_genes") {
      # getting markers => calculating contamination factor
      if (method == "soupx:background_genes")
        markers = get_markers(marker_file, sample_id, 3)
      else
        markers = get_top_n_markers(marker_file, sample_id, sc, 25)
      
      useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = markers)
      sc = calculateContaminationFraction(sc,markers,useToEst=useToEst, forceAccept = T)
    }
    
    decont_matrix = adjustCounts(sc)
    cont_matrix = sc$toc
  } 
  
  # DECONTX
  else if (substring(method,0,7) == "decontx") {
    if (sample_id == "mouse_kidney") {
      filtered = read.csv(dir,header = TRUE,sep = "\t")
      cont_matrix = as.matrix(filtered)
    } else if (sample_id == "hgmm12k") {
      #cont_matrix = as.matrix(get_raw_hgmm(files$CellRanger, config$sample_ids)@assays$RNA@counts) 
      cont_matrix = as.matrix(get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts) 
    
      storage.mode(cont_matrix) <- "integer"
    }
    
    if (method == "decontx:with_cell_types") {
      ## Cell Annotations
      decont_matrix = decontX(cont_matrix, z=as.numeric(factor(cell_annotations)))$resList$estNativeCounts
    } else {
      decont_matrix = decontX(cont_matrix)$resList$estNativeCounts
    }
    
  } 
  # CELLBENDER
  else if (method == "cellbender") {
    if (files$special[i] == FALSE) {
      stop("No special_files given, but path to filtered files is required for no decontamination")
    }
    
    if (sample_id != "hgmm12k") {
      filtered = read.csv(files$Filtered[i],header = TRUE,sep = "\t")
      cont_matrix = as.matrix(filtered)
    } else if (sample_id == "hgmm12k") {
      cont_matrix = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts 
    }
    if (config$run_cellbender) {
      # TODO: Fix for mouse_kidney dataset (specifically files$CellRanger)
      input_dir = paste(head(str_split(dir,"/")[[1]],-1),collapse="/") # removes the file name
      cellbender_args = c("remove-background", "--input", files$CellRanger, "--output", input_dir,"--expected-cells", dim(cont_matrix)[2])
      system2("cellbender", cellbender_args)
    }
    
    decont_matrix <- Read10X_h5(dir,use.names=T)
    
    # formatting cell barcodes to be '<sample id>_<barcode>'
    dimnames(decont_matrix)[[2]] = sapply(dimnames(decont_matrix)[[2]], function(x) {paste(sample_id, substring(x, 0, nchar(x)-2), sep="_")}, USE.NAMES=F)
    
    # only keep filtered cells and genes
    decont_matrix = decont_matrix[,(dimnames(decont_matrix)[[2]] %in% dimnames(cont_matrix)[[2]])]
    decont_matrix = decont_matrix[(dimnames(decont_matrix)[[1]] %in% dimnames(cont_matrix)[[1]]),]
    
    # reformatting cell barcodes to match barcodes in cell_annotations
    l = sapply(str_split(colnames(decont_matrix),"_"), function(n) paste(tail(n,1),"-1",sep=""))
    colnames(decont_matrix) = l
  } 
  # NO DECONTAMINATION
  else {
    # not actually "decontaminated" - however named this way
    if (sample_id == "mouse_kidney") {
      filtered = read.csv(dir,header = TRUE,sep = "\t")
      decont_matrix = as.matrix(filtered)
    } else if (sample_id == "hgmm12k") {
      decont_matrix = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts 
    }    
  } 
  
  out$seurat = CreateSeuratObject(decont_matrix)
  
  if (method == "cellbender") {
    # "fixing" cell annotations
    cell_annotations = cell_annotations[names(cell_annotations) %in% names(Idents(out$seurat))]
    cell_annotations = cell_annotations[order(match(names(cell_annotations), names(Idents(out$seurat))))]
  }
  # Add cluster information
  if (method != "none")
    Idents(out$seurat) <- cell_annotations

  return(out)
}


################################################################################################
# Removes sample_id (if present) and adds -1 (if not present) - for consistency
################################################################################################
fix_barcodes <- function(seurat) {
  # removing sample id
  if (str_count(colnames(seurat)[1],"_") > 0) {
    z <- sapply(colnames(seurat), function(x) {
      barcode = paste(tail(strsplit(x, "_")[[1]],1),collapse="_")
      return(barcode)
    }, USE.NAMES = F)
    
    log_print(paste("Sample ID being removed from barcodes\n", 
                      "From: ", colnames(seurat)[1],
                      "\nTo: ", z[1], sep=""))
    
    seurat = Seurat::RenameCells(seurat, new.names=z)
  }
  
  # adding -1
  if (str_count(colnames(seurat)[1],"-") == 0) {
    z <- sapply(colnames(seurat), function(x) {
      barcode = paste(x, "-1", sep="")
      return(barcode)
    }, USE.NAMES = F)
    
    log_print(paste("'-1' being added to barcodes\n", 
                      "From: ", colnames(seurat)[1],
                      "\nTo: ", z[1], sep=""))
    
    seurat = Seurat::RenameCells(seurat, new.names=z)
  }
  
  return(seurat)
}



################################################################################################
# Creates a tsv for each `sc.decont$toc` object in samples - prior to merging/integration 
################################################################################################
save_matrices <- function(samples) {
  file_dir = paste(files$output, "/matrices/",sep="")

  lapply(samples, function(x) {
    print(paste("Saving",x$sample_id))
	  write.table(as.matrix(x$seurat@assays$RNA@counts), file=paste(file_dir,x$sample_id,".tsv",sep=""),quote=FALSE,sep="\t")
  })
}



################################################################################################
# Returns a list of each cell type and the marker genes associated to each cell type
################################################################################################
get_markers <- function(dir_or_seurat, dataset, n=FALSE, b_compress=FALSE) {
  # n is the number of marker genes per cell type to select
  if (dataset != "hgmm12k") {
    suppressMessages({
      # order of celltypes in excel file
      excel_order = c("aLOH","B","CD_IC","CD_PC","CD_Trans","DCT-CNT","Endo","Fib","MC","MPH","NK","Podo","PT","T")
      markers = list()
      all_markers = c()

      # looping through sheets in excel file
      for (i in seq(length(excel_order))) {
        x = c(read_excel(dir_or_seurat, sheet=i,col_names=FALSE))$...1
        attr(x, "names") <- NULL

        # don't compress
        if (!b_compress) {
          if (n == FALSE) {
            markers[[i]] <- x
          } else {
            markers[[i]] <- x[1:n]
          }
        } 
        
        # compress into single vector w/ names being ct
        else {
          names(x) <- rep(excel_order[i],length(x))
          all_markers = c(all_markers, x)
        }
      }

      
      if (!b_compress)
        names(markers) = excel_order
    })
  } else if (dataset == "hgmm12k") {
    # Find markers
    raw_markers = FindMarkers(dir_or_seurat, ident.1 = "hg19", ident.2 = "mm10", verbose = F, logfc.threshold = 1, min.pct=0.5,test.use = "wilcox")
    
    # filter markers
    raw_markers = raw_markers[raw_markers$p_val_adj < 0.05,]
    
    # produce output
    markers = list(
      "hg19" = rownames(raw_markers[raw_markers$avg_logFC > 0,]),
      "mm10" = rownames(raw_markers[raw_markers$avg_logFC < 0,])
    )
    
    # compress into single vector w/ names being ct
    if (b_compress) {
      all_markers = c(markers$hg19, markers$mm10)
      names(all_markers) = c(rep("hg19", length(markers$hg19)), rep("mm10", length(markers$mm10)))
    }
  }
  
  return(if (b_compress) all_markers else markers)
}



################################################################################################
# Filters all the given marker genes (from `get_markers`) by selecting the top `n` in regards to `sc$soupProfile$est`
################################################################################################
get_top_n_markers <- function(dir, dataset, sc, n) {
  # getting the marker genes from the excel file 
  markers = get_markers(dir, dataset)
  all_markers = get_markers(dir, dataset, b_compress=TRUE) # a single vector containing all genes
  
  all_markers.no_dups <- all_markers[!all_markers %in% names(which(table(all_markers)!=1))]#all_markers[which(table(all_markers)==1)] #removing markers for >1 cell
  
  # creating `genes_ordered` -> [gene: `name`, rank: `numerical rank`]
  soupprofile = sc$soupProfile[rownames(sc$soupProfile) %in% all_markers.no_dups,]
  genes_ordered =soupprofile[order(soupprofile$est,decreasing=TRUE),]
  genes_ordered = data.frame("gene"=(rownames(genes_ordered)), "rank"=1:length(genes_ordered$est))
  genes_ordered$gene = unfactor(genes_ordered$gene)
  
  # looping through markers - filtering out 
  markers_top <- lapply(markers, function(x) {
    x <- lapply(x, function(y) {
      if (length(genes_ordered$rank[genes_ordered$gene==y]) != 0) {
        if (genes_ordered$rank[genes_ordered$gene==y]<=n) {
          return(y)
        }
      }
    })
    x <- unlist(x[lengths(x)!=0])
  })
  
  markers_top = markers_top[lengths(markers_top) != 0]
  return(markers_top)
}
