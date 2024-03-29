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
  
  ####################
  # SOUPX
  ####################
  if (substring(method,0,5)=="soupx") {
    # loads dir in 'SoupChannel' object
    if (sample_id != "hgmm12k") {
      sc = load10X(dir)

      marker_file = files$special[i] # setting marker file as dir to gene_sig file
    }

    else if (sample_id == "hgmm12k" & !config$stability_testing) {
      filtered = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)
      marker_file = filtered # setting marker file as seurat object
      
      # converting seurat objects to sparse matrices to create SC obj
      filtered = filtered@assays$RNA@counts
      raw = get_raw_hgmm(files$CellRanger, config$sample_ids)@assays$RNA@counts
      sc = SoupChannel(raw, filtered)
    }
	
    # adding clusters to sc obj
    sc = setClusters(sc, setNames(cell_annotations, names(cell_annotations))) 
    
    # removing any genes in the SC that are not in the filtered data
    if (length(sc$soupProfile[,1]) != length(rownames(sc$toc)))
      sc$soupProfile = sc$soupProfile[rownames(sc$toc),]
    
    # AUTOESTCONT
    if (method == "soupx:autoEstCont") {
      if (config$soupx_auto_vary_params == FALSE)
        sc = autoEstCont(sc)
      else {
        print(paste("Trying", config$soupx_auto_vary_params, 0.75))
        sc = try_autoEstCont(sc, 0.75, 0.25)
      }

    # MARKERS
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



  ####################
  # FastCAR
  ####################
  else if (method=="fastcar") {
    if (sample_id != "hgmm12k") {
      # needs gzipped files (ie barcodes.tsv.gz, features.tsv.gz, etc.)
      cell = read.cell.matrix(paste(dir, "filtered_gene_bc_matrices/mm10", sep="/")) # filtered
      full = read.full.matrix(paste(dir, "raw_gene_bc_matrices/mm10", sep="/")) # raw
    }

    else if (sample_id == "hgmm12k" & !config$stability_testing) {
      cell = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts # filtered
      full = get_raw_hgmm(files$CellRanger, config$sample_ids)@assays$RNA@counts # raw
    }
	
    # mm kidney emptyDropletCutoff, contaminationChanceCutoff = 100,0.05
    ambientProfile = determine.background.to.remove(full, cell, 300, 0.03)
    decont_matrix  = remove.background(cell, ambientProfile)
  } 



  ####################
  # DECONTX
  ####################
  else if (substring(method,0,7) == "decontx") {
    print(config$dataset)
	  if (config$dataset == "mouse_kidney") {
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
    } else if (method == "decontx:no_cell_types") {
      decont_matrix = decontX(cont_matrix)$resList$estNativeCounts
    } else if (method == "decontx:paper") {
      decont_matrix = decontX(cont_matrix, z=as.numeric(factor(cell_annotations)), maxIter = 60)$resList$estNativeCounts
    }
  } 

  ####################
  # CELLBENDER
  ####################
  else if (method == "cellbender") {
    if (files$special[i] == FALSE) 
      stop("No special_files given, but path to filtered files is required for no decontamination")
    
    # Data input
    if (sample_id != "hgmm12k") {
      filtered = read.csv(files$Filtered[i],header = TRUE,sep = "\t")
      cont_matrix = as.matrix(filtered)
    } else if (sample_id == "hgmm12k") {
      cont_matrix = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts 
    }

    # whether to run cellbender OR just read in results
    if (config$run_cellbender) {
      # TODO: check for mouse kidney data
      if (sample_id != "hgmm12k") 
        input_dir  = paste(files$CellRanger, sample_id, "raw_gene_bc_matrices", "mm10", sep="/")
      else if (sample_id == "hgmm12k")
        input_dir  = files$CellRangerMerged

      output_dir = paste(files$output, sample_id, sep="/")
      filename   = paste(sample_id, ".h5", sep="")
      
      # checks for output dir presense - if it doesnt exist, cellbender wont realise until after its finished executing
      if (!dir.exists(output_dir))
        dir.create(output_dir)
      
      cellbender_args = c("remove-background", "--input", input_dir, 
                          "--output", paste(output_dir, filename, sep="/"),
                          "--expected-cells", dim(cont_matrix)[2],
      			  if (config$cellbender_gpu) "--cuda" else NULL)
      
      print("##########")
      print("The arguements for CellBender are:")
      print(cellbender_args)
      print("##########")

      system2("cellbender", cellbender_args)
    }

    if (config$benchmarking)
      decont_matrix <- Read10X_h5(paste(output_dir, "/", sample_id, "_filtered.h5", sep=""), use.names=T)
    else
      decont_matrix <- Read10X_h5(paste(paste(files$output, sample_id, sample_id, sep="/"), "_filtered.h5", sep=""), use.names=T)
    
    # formatting cell barcodes to be '<sample id>_<barcode>'
    if (sample_id != "hgmm12k")
      dimnames(decont_matrix)[[2]] = sapply(dimnames(decont_matrix)[[2]], function(x) {paste(sample_id, substring(x, 0, nchar(x)-2), sep="_")}, USE.NAMES=F)
    else # converts to seurat object then back to matrix to convert `_` to `-`
      decont_matrix = CreateSeuratObject(decont_matrix)@assays$RNA@counts
    
    # only keep filtered cells and genes
    decont_matrix = decont_matrix[,(dimnames(decont_matrix)[[2]] %in% dimnames(cont_matrix)[[2]])]
    decont_matrix = decont_matrix[(dimnames(decont_matrix)[[1]] %in% dimnames(cont_matrix)[[1]]),]
    
    # reformatting cell barcodes to match barcodes in cell_annotations
    if (sample_id != "hgmm12k")
      colnames(decont_matrix) = sapply(str_split(colnames(decont_matrix),"_"), function(n) paste(tail(n,1),"-1",sep=""))
  } 

  ####################
  # NO DECONTAMINATION
  ####################
  else {
    # not actually "decontaminated" - however named this way
    if (sample_id != "hgmm12k") {
      decont_matrix = Read10X(paste(dir, "filtered_gene_bc_matrices", "mm10", sep="/"))#@assays$RNA@counts
      
      # old way:
      #filtered = read.csv(special_files, header = TRUE, sep = "\t")
      #decont_matrix = as.matrix(filtered)
    } else if (sample_id == "hgmm12k") {
      decont_matrix = get_filtered_hgmm(files$CellRanger, files$CellAnnotations, config$sample_ids)@assays$RNA@counts 
    }    
  } 
  
  ####################
  # FINAL PROCESSING
  ####################
  # create seurat obj
  out$seurat = CreateSeuratObject(decont_matrix)
  
  # "fixing" cell annotations for cell bender
  if (method == "cellbender") {
    cell_annotations = cell_annotations[names(cell_annotations) %in% names(Idents(out$seurat))]
    cell_annotations = cell_annotations[order(match(names(cell_annotations), names(Idents(out$seurat))))]
  }

  # Add cluster information to seurat object
  Idents(out$seurat) <- cell_annotations

  return(out)
}


################################################################################################
# Recursive function to test soupx:autoestcont params until it works
################################################################################################
try_autoEstCont <- function(sc, val, iter) {
  return(tryCatch({
    # check which param(s) to vary / use
    # try to calculate with soupx
    if (config$soupx_auto_vary_params == "tfidfMin")
      sc = autoEstCont(sc, tfidfMin = (val))
    else if (config$soupx_auto_vary_params == "soupQuantile")
      sc = autoEstCont(sc, soupQuantile = (val))
    else
      sc = autoEstCont(sc, tfidfMin = val[1], soupQuantile = val[2])

    # on success - print message and return soupchannel object
    if (config$soupx_auto_vary_params == "tfidfMin" || config$soupx_auto_vary_params == "soupQuantile")
      print(paste("soupx:autoEstCont successful w/ ", config$soupx_auto_vary_params, "=", val))
    else 
      print(paste("soupx:autoEstCont successful w/ tfidfMin and soupQuantile =", val))

    return(sc)
  }, error = function(e) { # on error - print message & recall this function w/ `val-iter`
    if (config$soupx_auto_vary_params == "tfidfMin" || config$soupx_auto_vary_params == "soupQuantile")
      if ((val-iter) >= 0) # check val-iter >= 0
        print(paste("Reducing", config$soupx_auto_vary_params, "to", val-iter))
      else 
        stop(paste("varying", config$soupx_auto_vary_params, "did not work (still erroring @", val,")"))
    else 
      if ((val-iter) >= 0) # check val-iter >= 0
        print(paste("Reducing tfidfMin and soupQuantile to", val-iter))
      else 
        stop(paste("varying tfidfMin and soupQuantile did not work (still erroring @ 0,0)"))

    try_autoEstCont(sc, val-iter, iter)
  }))
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
  if (config$benchmarking) { # different directories for benchmarking and normal operation
    file_dir = paste(files$output, sample_id, "matrices/", sep="/")
    if (!dir.exists(file_dir))
      dir.create(file_dir)
  } else
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
    raw_markers = FindMarkers(dir_or_seurat, ident.1 = "hg19", ident.2 = "mm10", verbose = F, logfc.threshold = 1, min.pct=0.5, test.use = "wilcox")
    
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
}
