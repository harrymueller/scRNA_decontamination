
################################################################################################
# Re-clusters the given seurat object using gene signatures & Seurat:AddModuleScore
# Uses the Monte-Carlo procedure and adjusts the p-value via Benjamini & Hochberg 
################################################################################################
reCluster <- function(seurat) {
  gene_sig_file = files$GeneSignatures
  alpha = config$alpha

  metadata_length_pre = length(names(seurat@meta.data))
  gene_sigs = get_markers(gene_sig_file, config$dataset)
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
# Returns a character vector of cell annotations
    # if use_new -> reads from a tsv of the same format as reclustered cell annotations are saved as
    # else -> reads a xlsx file of the same format as the cell annotations from the kidney paper
##########################################
get_clusters <- function(path, sample_id, is_xlsx) {
  if (sample_id == "hgmm12k") {
    classifications = read.csv(path)
    classifications = classifications[classifications$call != "Multiplet",] # removing multiplets
    
    cell_annotations = c(unfactor(classifications$call))
    names(cell_annotations) = classifications$barcode
    return(cell_annotations)
  }
  
  else {
    if (is_xlsx) {
      is_preserved = length(str_split(sample_id, "_",simplify=TRUE))>2

      # Annotations for preserved cells are on a separate sheet within the file
      if (!is_preserved) {
        suppressMessages({
          cell_annotations = read_excel(path, sheet=1)
        })
        # library names are different from preserved library names
        cell_annotations = cell_annotations[cell_annotations$Library==sample_id,]

      } else {
        suppressMessages({
          cell_annotations = read_excel(path, sheet=2)
        })
        cell_annotations = cell_annotations[cell_annotations$Preservation == "MeOH",c(1:3,5)]

        cell_annotations = cell_annotations[cell_annotations$Library==str_split(sample_id, "_",simplify=TRUE)[1],]
      }

      names(cell_annotations)[1] = "barcodes"
      names(cell_annotations)[4] = "orig.ident"

      cell_annotations$barcodes = sapply(cell_annotations$barcodes, function(x) paste(substring(x, nchar(sample_id)+2), "-1", sep=""))
      names(cell_annotations$orig.ident) = cell_annotations$barcodes

      #cell_annotations = subset(cell_annotations, select=-c(1,2,3))
      return(cell_annotations$orig.ident)
    } else {
      annotations = read.table(path, sep="\t",header=F,skip=1,stringsAsFactors=F, as.is=T)

      names(annotations) = c("barcode", "celltype")

      # creating sample_id variable
      annotations$sample_id = sapply(annotations$barcode, FUN=function(x) {
        return(paste(head(str_split(x,"_",)[[1]],-1),collapse = "_"))
      })

      annotations = annotations[which(annotations$sample_id==sample_id),]


      # removing sample_id from barcode
      annotations$barcode = sapply(annotations$barcode, FUN=function(x) {
        return(tail(str_split(x, "_")[[1]],1))
      })

      ct = annotations$celltype
      names(ct) = annotations$barcode
      return(ct)
    }
  }
}
