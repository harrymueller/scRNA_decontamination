################################################################################################
## Gets all the DEGs and returns a list object
################################################################################################
get_DEGs <- function(samples.combined, save_name=FALSE) {
  c1 = (unique(samples.combined@meta.data$celltype[which(samples.combined@meta.data$preservation=="fresh")]))
  c2 = (unique(samples.combined@meta.data$celltype[which(samples.combined@meta.data$preservation=="MeOH")]))

  cell_types = c1[which(c1%in% c2)]
  cell_types = cell_types[which(cell_types != "Unknown")] #removing unknown cell type

  DEGs = as.list(cell_types) 
  names(DEGs) <- cell_types

  for (ct in cell_types) {
    # logfc.threshold = 1
    tryCatch({
      DEGs[[ct]] <- FindMarkers(samples.combined, ident.1 = paste(ct,"fresh",sep="_"), ident.2 = paste(ct,"MeOH",sep="_"), 
                                verbose = F, logfc.threshold = 1, min.pct=0.5,test.use = "wilcox")
    }, error = function(e) {
      print(e)
      log_print(paste(e))
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
	
  if (save_name != FALSE) 
    saveRDS(DEGs,save_name)

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
DEGs_histogram <- function(DEGs) {
  f_name = paste(files$output,"/plots/DEG_histograms.png",sep="")
  ct_order = config$ct_order_dotplots

  DEGs.c = do.call("rbind", DEGs)
  ct_order = ct_order[which(ct_order %in% names(DEGs))]
  
  DEGs.c$cell_type <- ordered(DEGs.c$cell_type,ct_order)
  DEGs.c = DEGs.c[order(DEGs.c$cell_type, decreasing=F),]

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

  if (config$fonts) {
    p1 = p1 + theme(text=element_text(size=16, family="TT Times New Roman"))
    p2 = p2 + theme(text=element_text(size=16, family="TT Times New Roman"))
  }

  p <- p1 + p2
                    
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
  
  # order of ct for dotplot
  order_paper = order_paper[which(order_paper %in% levels(samples.combined))]
  s=20 # text size

  # ct not included in dot plot
  not_inc = unique(Idents(samples.combined))[which(!(unique(Idents(samples.combined)) %in% order_paper))]

  # Expression in FRESH
  # celltype.fresh <- normal ct for fresh sample | ct_meoh for rest <- selects only fresh
  Idents(samples.combined) <- "celltype.fresh"

  # changing levels of samples.combined to be correct order
  l = c(order_paper, paste(order_paper,"MeOH",sep="_"), paste(not_inc, "", sep=""), paste(not_inc, "MeOH", sep="_"), use.names=F)
  
  levels(samples.combined) <- l
  # genes from paper dotplot
  p1 <- DotPlot(samples.combined,idents=order_paper,features = genes_paper, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s, angle = 60, hjust = 1)) + ggtitle("A. DEGs detected in the paper;\n   Expression in fresh samples after decontamination")
  # DEGs from this analysis dotplot
  p3 <- DotPlot(samples.combined,idents=order_paper,features = c(genes.under,genes.over), cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("B. DEGs detected after decontamination;\n   Expression in fresh samples after decontamination")
	
  # Expression in MEOH
  Idents(samples.combined) <- "celltype.meoh"

  # changing levels of samples.combined to be correct order
  l = c(order_paper, paste(order_paper,"fresh",sep="_"), paste(not_inc, "", sep=""), paste(not_inc, "fresh", sep="_"), use.names=F)
  levels(samples.combined) <- l

  # genes from paper dotplot
  p2 <- DotPlot(samples.combined,idents=order_paper,features = genes_paper, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("C. DEGs detected in the paper;\n   Expression in MeOH samples after decontamination")
  # degs from this analysis dotplot
  p4 <- DotPlot(samples.combined,idents=order_paper,features = c(genes.under, genes.over), cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + xlab("DEGs") + ylab("Cell Type") + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + ggtitle("D. DEGs detected after decontamination;\n   Expression in MeOH samples after decontamination")
  
  p <- (p1+p2)/(p3+p4)

  ggsave(f_name,p,width=16, height=12)
}

# separate dot plots for under and over expr in MeOH
# about 5-10 -> include num of ct req.
# + under expr >= 9ct
DEGs_dotplot_specific <- function(samples.combined, order_paper, output_dir, DEGs.over, DEGs.under) {
  Idents(samples.combined) <- "celltype"
  s=20 # text size
  num_ct = length(unique(Idents(samples.combined)))

  # order of ct for dotplot
  order_paper = order_paper[which(order_paper %in% levels(samples.combined))]

  # plot DEGs over expr in meoh
  genes.over = degs_opt_ct(DEGs.over, num_ct)
  genes.over.n = genes.over[1]
  genes.over = genes.over[seq(2, length(genes.over))]

  # celltype.fresh <- normal ct for fresh sample | ct_meoh for rest <- selects only fresh
  Idents(samples.combined) <- "celltype.fresh" 
  p1 = DotPlot(samples.combined, idents=order_paper, features = genes.over, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("DEGs Over-Expressed in MeOH\n>= ", genes.over.n, " Celltypes\nFresh Cells"))
  Idents(samples.combined) <- "celltype.meoh"
  p2 = DotPlot(samples.combined, idents=order_paper, features = genes.over, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("MeOH Cells"))
  ggsave(paste(output_dir,"/DEG_higher_expression_MeOH.png",sep=""),
    p1/p2,width=10, height=14)

  # plot DEGs under expr in MeOH
  genes.under = degs_opt_ct(DEGs.under, num_ct)
  genes.under.n = genes.under[1]
  genes.under = genes.under[seq(2, length(genes.under))]

  Idents(samples.combined) <- "celltype.fresh" 
  p1 = DotPlot(samples.combined, idents=order_paper, features = genes.under, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("DEGs Under-Expressed in MeOH\n>= ", genes.under.n, " Celltypes\nFresh Cells"))
  Idents(samples.combined) <- "celltype.meoh" 
  p1 = DotPlot(samples.combined, idents=order_paper, features = genes.under, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("MeOH Cells"))
  ggsave(paste(output_dir,"/DEG_lower_expression_MeOH.png",sep=""),
    p1/p2,width=10, height=14)

  # under expr in meoh in >= 9ct
  genes.nine.under = get_genes_de(DEGs.under, 8) # >= 9

  Idents(samples.combined) <- "celltype.fresh" 
  p1 = DotPlot(samples.combined, idents=order_paper, features = genes.nine.under, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("DEGs Under-Expressed in MeOH\n>= ", 9, " Celltypes\nFresh Cells"))
  Idents(samples.combined) <- "celltype.meoh" 
  p2 = DotPlot(samples.combined, idents=order_paper, features = genes.nine.under, cols="RdYlBu", dot.scale = 5,cluster.idents = F, scale.min=0, scale.max=100) + 
    xlab("DEGs") + ylab("Cell Type") + 
    theme(text = element_text(size=s+2), axis.text.y = element_text(size=s), axis.text.x = element_text(size=s,angle = 60, hjust = 1)) + 
    ggtitle(paste0("MeOH Cells"))
  ggsave(paste(output_dir,"/DEG_lower_expression_9_MeOH.png",sep=""),
    p1/p2,width=10, height=14)
  
  Idents(samples.combined) <- "celltype"
}

# loop through number of ct to find num degs between 7 and 10
degs_opt_ct <- function(DEGs, num_ct) {
  for (i in seq(num_ct)) {
    genes <- get_genes_de(DEGs, i)
    if (length(genes) > 7 && length(genes) <= 10) return(c(i, genes))
  }
  return(c(2, get_genes_de(DEGs, 2)))
}