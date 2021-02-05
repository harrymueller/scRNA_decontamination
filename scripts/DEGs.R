################################################################################################
### Adds some metadata to the combined seurat object
################################################################################################
adding_metadata <- function(samples.combined) {
  #order_paper <- rev(c("Fib", "MPH", "Podo", "aLOH", "Unknown","CD_PC", "T","NK","DCT","CNT","B","MC","CD_IC","Endo","PT","CD_Trans"))
  Idents(samples.combined)[which(Idents(samples.combined)=="CD_Trans")] = "Unknown"
  order_paper = order_paper[which(order_paper %in% levels(samples.combined))]
  
  #samples.combined <- readRDS(paste(output, "samples_integrated_rd.Rda", sep="/"))
  DefaultAssay(samples.combined) <- "RNA"
  
  # adding metadata for `celltype_method`, `celltype`, and changing default ident to `celltype_method`
  samples.combined$celltype_method <- paste(Idents(samples.combined), samples.combined$preservation, sep = "_")
  samples.combined$celltype <- Idents(samples.combined)
  Idents(samples.combined) <- "celltype_method"
  
  samples.combined$celltype.fresh = unlist(lapply(samples.combined$celltype_method, function(x) {
    end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="fresh")
    x <- if (end) str_sub(x,end=-7) else x
  }))
  samples.combined$celltype.meoh = unlist(lapply(samples.combined$celltype_method, function(x) {
    end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="MeOH")
    x <- if (end) str_sub(x,end=-6) else x
  }))
  
  return(samples.combined)
}



################################################################################################
## Gets all the DEGs and returns a list object
################################################################################################
get_DEGs <- function(samples.combined, save_name=FALSE) {
  #cell_types = varhandle::unfactor(samples.combined$celltype)
  #cell_types = names(which(table(cell_types)>10)) #filtering out cell types with less than 10 cells
  c1 = (unique(samples.combined@meta.data$celltype[which(samples.combined@meta.data$preservation=="fresh")]))
  c2 = (unique(samples.combined@meta.data$celltype[which(samples.combined@meta.data$preservation=="MeOH")]))
  cell_types = c1[which(c1%in% c2)]
  cell_types = cell_types[which(cell_types != "Unknown")] #removing unknown cell type
  DEGs = as.list(cell_types)
  names(DEGs) <- cell_types

  for (ct in cell_types) {
    # logfc.threshold = 1
    tryCatch({
	  # TODO catch < 3 cells prior to error > error now is no DEGs
	  DEGs[[ct]] <- FindMarkers(samples.combined, ident.1 = paste(ct,"fresh",sep="_"), ident.2 = paste(ct,"MeOH",sep="_"), 
                              verbose = F, logfc.threshold = 1, min.pct=0.5,test.use = "wilcox")
	}, error = function(e) {
		print(e)
	  log_warning(paste(e), paste(output, "errors.txt", sep="/")) #TODO improve
      DEGs[ct] = 0
	}) 
  }

  # adding cell type to individual df and filtering out genes w/ FDR > 0.05
  DEGs = lapply(names(DEGs), function(x) {
	if (is.null(dim(DEGs[[x]]))) {
		DEGs[[x]] <- NULL
	} else {
	    adjusted = DEGs[[x]][DEGs[[x]]$p_val_adj<0.05,]

	    if (length(adjusted[,1]) > 0) {
	      adjusted$cell_type = factor(x)
	    } 

		adjusted$gene = rownames(adjusted)
	    DEGs[[x]] = adjusted
	}
  })
  names(DEGs) <- cell_types
	
  if (save_name != FALSE) {
    saveRDS(DEGs,save_name)
  }
  #DEGs = readRDS(paste(output, "DEGs_per_celltype.Rda", sep="/"))
  return(DEGs)
}



################################################################################################
### DEGs over/underexpressed in MeOH cells
################################################################################################
get_over_under_DEGs <- function(DEGs, b_over) {
  # filtering only genes over/underexpressed in MeOH samples
  DEGs.selected = DEGs

  DEGs.selected = lapply(names(DEGs.selected), function(x) {
    if (b_over) {
      DEGs.selected[[x]] <- DEGs.selected[[x]][DEGs.selected[[x]]$avg_logFC < 0,]
    } else {
      DEGs.selected[[x]] <- DEGs.selected[[x]][DEGs.selected[[x]]$avg_logFC > 0,]
    }
  })
  
  names(DEGs.selected) <- names(DEGs)
  return(DEGs.selected)
}



################################################################################################
#### Number of DEGs in at least `num_cells` cell types
# returns the names of genes DE in at least `num_cells`
################################################################################################
get_genes_de <- function(DEGs, num_cells) {
  DEGs.c = do.call("rbind", DEGs)

  # Genes DE in at least 8 cells
  genes = names(which(table(DEGs.c$gene)>num_cells))
  return(genes)
}



################################################################################################
#### Outputting DEGs to `xlsx` document
################################################################################################
save_DEGs <- function(DEGs, f_name, genes.over, genes.under) {
  wb <- createWorkbook()
  
  # DEGs over/under expressed >9 cell types
  s <- createSheet(wb, "overexpressed_genes_>8_celltypes")
  addDataFrame(c(
    if (is.null(genes.over)) NULL else c("Genes over-expressed in MeOH cells in at least 9 cells types",genes.over),
    if (is.null(genes.under)) NULL else c("Genes under-expressed in MeOH cells in at least 9 cells types",genes.under)
  ), s, row.names = FALSE, col.names = FALSE)
  
  # DEGs by cell type
  lapply(names(DEGs), function(x) {
    s <- createSheet(wb, x)
    addDataFrame(DEGs[[x]], s)
  })
  
  saveWorkbook(wb, f_name)
}



################################################################################################
### barplots of number of DEGs by (celltype, over/under expression)
################################################################################################
DEGs_histogram <- function(DEGs, f_name) {
  DEGs.c = do.call("rbind", DEGs)
  ct_order = rev(c("Fib", "MPH", "Podo", "aLOH", "Unknown","CD_PC", "T","NK","DCT","CNT","DCT_CNT","B","MC","CD_IC","Endo","PT"))
  ct_order = ct_order[which(ct_order %in% names(DEGs))]
  
  DEGs.c$cell_type <- ordered(DEGs.c$cell_type,ct_order)
  
  DEGs.c = DEGs.c[order(DEGs.c$cell_type, decreasing=F),]
  #DEGs.c$gene = factor(sapply(rownames(DEGs.c), function(x) str_split(x,fixed("."),simplify = T)[2]))
  
  #temp = data.frame(table(DEGs.c$cell_type), DEGs.c$)
  #names(temp) <- c("celltype", "degs")
  p1 <-ggplot(data.frame(table(DEGs.c[DEGs.c$avg_logFC < 0,]$cell_type)), aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") +
    coord_flip() + ylab("Number of DEGs") +
    scale_x_discrete(name="Cell Type") + ggtitle("DEGs up-regulated in MeOH samples") + 
    (if (max(sapply(DEGs,function(x) length(x$avg_logFC[x$avg_logFC<0]),USE.NAMES = F)) < 21) ggplot2:::limits(c(0,20),"y"))
  
  p2<- ggplot(data.frame(table(DEGs.c[DEGs.c$avg_logFC > 0,]$cell_type)), aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") +
    coord_flip() + ylab("Number of DEGs") +
    scale_x_discrete(name="Cell Type") + ggtitle("DEGs down-regulated in MeOH samples") + 
    (if (max(sapply(DEGs,function(x) length(x$avg_logFC[x$avg_logFC>0]),USE.NAMES = F)) < 21) ggplot2:::limits(c(0,20),"y"))
  p <- p1 + p2 #+ plot_annotation(
  #  title="Differentially Expressed Genes from Methanol-Fixed Samples Compared to Fresh Samples"
  #)
#  p
  ggsave(f_name,p,width=9, height=5)
}


num_DEGs <- function() {
  # unique
  length(get_genes_de(DEGs.over,0))
  length(get_genes_de(DEGs.under,0))
  # total
  length(rownames(do.call("rbind",DEGs.over)))
  length(rownames(do.call("rbind",DEGs.under)))
}



################################################################################################
##### overexpressed >= 9 cell types <- dotplots
################################################################################################
DEGs_dotplot_over_under_expression <- function(samples.combined, f_name, order_paper, genes_paper, genes.over, genes.under) {
  Idents(samples.combined) <- "celltype"
  #Idents(samples.combined) = factor(Idents(samples.combined))
  order_paper = order_paper[which(order_paper %in% levels(samples.combined))]

  s=20
  Idents(samples.combined) <- "celltype.fresh"
  levels(samples.combined) <-  c(order_paper,paste(order_paper,"MeOH",sep="_"))

  p1 <- DotPlot(samples.combined,idents=order_paper,features = genes_paper, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s, angle = 60, hjust = 1)) + ggtitle("A. DEGs detected in the paper;\n   Expression in fresh samples after decontamination")
  p3 <- DotPlot(samples.combined,idents=order_paper,features = c(genes.under,genes.over), cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("B. DEGs detected after decontamination;\n   Expression in fresh samples after decontamination")
  
  Idents(samples.combined) <- "celltype.meoh"
  levels(samples.combined) <-  c(order_paper,paste(order_paper,"fresh",sep="_"))
  p2 <- DotPlot(samples.combined,idents=order_paper,features = genes_paper, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("C. DEGs detected in the paper;\n   Expression in MeOH samples after decontamination")
  p4 <- DotPlot(samples.combined,idents=order_paper,features = c(genes.under, genes.over), cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("D. DEGs detected after decontamination;\n   Expression in MeOH samples after decontamination")
  
  p <- (p1+p2)/(p3+p4)

  ggsave(f_name,p,width=16, height=12)
}
