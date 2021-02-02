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
  } 
  else if (substring(method,0,7) == "decontx") {
    filtered = read.csv(dir,header = TRUE,sep = "\t")
    cont_matrix = as.matrix(filtered)
    
    if (method == "decontx:with_cell_types") {
      ## Cell Annotations
      cell_annotations = as.numeric(factor(cell_annotations))
      
      decont_matrix = decontX(cont_matrix, z=cell_annotations)$resList$estNativeCounts
    } else {
      decont_matrix = decontX(cont_matrix)$resList$estNativeCounts
    }
    
  } 
  else if (method == "cellbender") {
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
  } 
  else {
    filtered = read.csv(special_files, header = TRUE, sep = "\t")
    decont_matrix = as.matrix(filtered) # not actually "decontaminated" - however named this way
  } 
  
  out$seurat = CreateSeuratObject(decont_matrix)
  
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
save_matrices <- function(samples, file_dir) {
  # file_dir = "../data/output/no_decont/"
  lapply(samples, function(x) {
    print(paste("Saving",x$sample_id))
    write.table(as.matrix(x$seurat@assays$RNA@counts), file=paste(file_dir,x$sample_id,".tsv",sep=""),quote=FALSE,sep="\t")
  })
}


################################################################################################
# Re-clusters the given seurat object using gene signatures & Seurat:AddModuleScore
# Uses the Monte-Carlo procedure and adjusts the p-value via Benjamini & Hochberg 
################################################################################################
reCluster <- function(seurat, gene_sig_file, alpha) {
  # gene_sig_file = "/home/harry/gdrive/Education/Uni/GCRL2000/data/cell_type_gene_signatures.xlsx"
  metadata_length_pre = length(names(seurat@meta.data))
  gene_sigs = get_markers(gene_sig_file)
  names(gene_sigs)[which(names(gene_sigs) == "DCT-CNT")] = "DCT_CNT"
  
  # adding module scores per gene
  for (i in names(gene_sigs)) {
    seurat=AddModuleScore(seurat,gene_sigs[i],name=i,nbin=20, ctrl.size=50)
  }
  # renaming meta.data to match cts
  names(seurat@meta.data) = sapply(X=names(seurat@meta.data), FUN= function(x) {
    if (str_sub(x,-1,-1) == "1") {
      x <- str_sub(x, 0, -2)
    } else {
      x <- x
    }
  }, USE.NAMES=FALSE)
  
  cts = names(gene_sigs)
  
  # saves the two largest clusters and values
  clusters = c()
  values = c()
  clusters_2nd = c()
  values_2nd = c()
  for (i in seq(1:length(colnames(seurat)))) {
    clusters = append(clusters, cts[which.max(seurat@meta.data[i,cts])[[1]]])
    values = append(values, seurat@meta.data[i,clusters[i]])
    #values = append(values, seurat@meta.data[i, which.max(seurat@meta.data[i,cts])[[1]]+metadata_length_pre])
    
    # adding second best cluster to meta data
    cts_i = cts[which(cts != clusters[i])]
    
    clusters_2nd = append(clusters_2nd, cts_i[which.max(seurat@meta.data[i,cts_i])[[1]]])
    values_2nd = append(values_2nd, seurat@meta.data[i,clusters_2nd[i]])
  }
  
  # adds the difference of the top two module scores to the metadata
  seurat@meta.data$module_score_diff = values - values_2nd
  seurat@meta.data$module_score_ct = clusters_2nd
  
  # neg module score = unknown
  clusters[which(values<0)] = "Unknown"
  clusters[which(clusters=="CD_Trans")] = "Unknown"
  
  # splits the seurat object by the new clusters
  Idents(seurat) = clusters
  new.seurat = list()
  split = SplitObject(seurat)
  
  # loops through all subsets
  for (subset in split) {
    ct = levels(Idents(subset))
    # Skips unknown ct and and cts that only have 1 barcode
    print(ct)
    if (ct != "Unknown" && length(colnames(subset)) > 1) {
      ### FOLLOWING CODE IS ADAPTED FROM ELENA
      # create random gene sets of the same length as the signature for that cell type
      random_gene_sets <- lapply(vector("list", 1000), function(x){sample(rownames(subset), length(gene_sigs[ct][[1]]))})
      
      # calculate scores for these 1000 random gene sets
      nbin = 20 # default bin size
      repeat { # repeats until successfully (only required for some subsets) - if failed - halves bin size until success
        do_break = T
        tryCatch({
          subset <- AddModuleScore(subset, features = random_gene_sets, name = "RandomRun",nbin=nbin, ctrl.size = 50)
        }, error = function(e) {
          nbin = nbin/2
          
          # if bin size does decrease - put it in the log
          message("WARNING: Decreasing number of bins")
          log_print(paste("reCluster error:\nHalving number of bins to ", nbin,"\n",
                            "Sample ID = ", as.character(unique(samples.combined[[1]]@meta.data$orig.ident)), "\n",
                            "CT = ", ct, sep=""))
          do_break = F
        }, finally = function(e) {
          do_break = T
        })
        
        if (do_break == T) {
          break
        }
      }     
      
      # select columns with scores for 1000 random gene sets
      score_columns <- grep("RandomRun", colnames(subset@meta.data))
      # cbind together the scores from the real gene set and random gene sets 
      randomscores <- subset@meta.data[, score_columns]
      interim <- cbind(real = values[which(clusters==ct)], randomscores)
      # calculate how many random runs gave a score value above the real score (in the first column of 'interim')
      res <- apply(interim, 1, function(x){sum(x[-1] >= x[1])})
      # calculate p-value from Monte-Carlo procedure as described by North et al (they recommend (r+1)/(n+1) instead of r/n)
      res_pval <- (res+1)/1001
      # do BH adjustment
      res_fdr <- p.adjust(res_pval, "BH")
      
      clus = rep(ct, length(res_fdr))
      clus[which(res_fdr>alpha)] = "Unknown"
    } else if (ct != "Unknown" && length(colnames(subset)) == 1) {
      # if a ct only has 1 sample, cant use AddModuleScore <- so it is added to unknown
      message("WARNING: Number of cells within a new cluster = 1")
      log_print(paste("Number of cells within a cluster = 1:\nAdding barcode to unknown cluster\n",
                        "Sample ID = ", as.character(unique(samples.combined[[1]]@meta.data$orig.ident)), "\n",
                        "CT = ", ct, sep=""))
      clus = "Unknown"
    } 
    
    Idents(subset) = clus
    new.seurat = append(new.seurat, subset)
  }
  # Re-merging seurat object after splitting and checking p-vals
  new.seurat = merge(new.seurat[[1]], new.seurat[2:length(new.seurat)],add.cell.ids=NULL)
  
  # removing excess meta.data
  new.seurat@meta.data = new.seurat@meta.data[,which(!(names(new.seurat@meta.data) %in% cts))]
  new.seurat@meta.data = new.seurat@meta.data[,which(substring(names(new.seurat@meta.data), 0, 9) != "RandomRun")]
  
  return(seurat)
}



################################################################################################
# Returns a list of each cell type and the marker genes associated to each cell type
################################################################################################
get_markers <- function(dir, n=FALSE, b_compress=FALSE) {
  # n is the number of marker genes per cell type to select
  suppressMessages({
    # order of celltypes in excel file
    excel_order = c("aLOH","B","CD_IC","CD_PC","CD_Trans","DCT-CNT","Endo","Fib","MC","MPH","NK","Podo","PT","T")
    markers = list()
    all_markers = c()
    
    for (i in seq(length(excel_order))) {
      x = c(read_excel(dir, sheet=i,col_names=FALSE))$...1
      attr(x, "names") <- NULL
      
      if (!b_compress) {
        if (n == FALSE) {
          markers[[i]] <- x
        } else {
          markers[[i]] <- x[1:n]
        }
      } 
      
      else {
        names(x) <- rep(excel_order[i],length(x))
        all_markers = c(all_markers, x)
      }
    }
    
    if (!b_compress) {
      names(markers) = excel_order
    }
  })
  
  return(if (b_compress) all_markers else markers)
}

################################################################################################
# Filters all the given marker genes (from `get_markers`) by selecting the top `n` in regards to `sc$soupProfile$est`
################################################################################################
get_top_n_markers <- function(dir, sc, n) {
  # getting the marker genes from the excel file 
  markers = get_markers(dir)
  all_markers = get_markers(dir, b_compress=TRUE) # a single vector containing all genes
  
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