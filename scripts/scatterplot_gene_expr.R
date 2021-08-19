gene_expr_scatter_plots = function(prior, post) {
  DefaultAssay(prior) = "RNA"
  
  # for each pres
  for (pres in c("fresh", "MeOH")) {
    # get average expr
    prior_expr = AverageExpression(prior[,prior$preservation == pres], assay = "RNA")$RNA
    post_expr = AverageExpression(post[,post$preservation == pres], assay = "RNA")$RNA
    diff = (post_expr - prior_expr)
    # reorder
    prior_expr = prior_expr[order(match(rownames(prior_expr), rownames(post_expr))),]
    
    for (ct in names(post_expr)) {
      new = data.frame('X' = log(prior_expr[[ct]] + 1), 
                      'Y' = log(post_expr[[ct]] + 1), 
                      "NAME" = rownames(prior_expr))
      new["bLABEL"] = new$NAME %in% rownames(diff[order(diff[[ct]], decreasing = T),])[1:10]

      p = ggplot(new, aes(x=X, y=Y, label=NAME)) +  
        geom_point(alpha = 0.5) + ggtitle(paste("Average gene expression plotted on log(x+1) scales for", pres, "&", ct)) +
        geom_text_repel(aes(label=ifelse(bLABEL>0.2,as.character(NAME),'')),max.overlaps=100) + 
        geom_abline(slope=1, intercept=0) + 
        ylab("Post-Decontamination") + xlab("Pre-Decontamination") 
      
      ggsave(paste(files$output, "/plots/gene_expression/", ct, "_", pres, ".png", sep=""), p, width=8, height=8)
    }
  }
}

run_degs_prior_post <- function(prior, post) {
  prior@meta.data["type"] = "prior"
  post@meta.data["type"] = "post"

  merged = merge(prior, post, add.cell.ids = c("prior", "post"))
  DefaultAssay(merged) = "RNA" # just triple confirming

  # identifier for DEGs
  merged@meta.data["celltype_type"] = paste(merged$celltype, merged$type, sep="_")

  merged = SplitObject(merged, split.by = "preservation")
  for (pres in c("fresh", "MeOH")) {
    degs = get_DEGs_prior_post(merged[[pres]])

    summary = get_DEGs_summary(degs)

    save_DEGs_prior_post(degs, summary, paste(config$output_dir, "/", current_method, "/", pres, "_DEGs_between_prior_post.xlsx", sep=""))
  }
}

#canablised from DEGs.R
get_DEGs_prior_post <- function(samples) {
  Idents(samples) = "celltype_type"
  cell_types = unique(samples$celltype)
  cell_types = cell_types[which(cell_types != "Unknown")] #removing unknown cell type

  DEGs = as.list(cell_types) 
  names(DEGs) <- cell_types

  for (ct in cell_types) {
    # logfc.threshold = 1
    tryCatch({
      DEGs[[ct]] <- FindMarkers(samples, ident.1 = paste(ct,"prior",sep="_"), ident.2 = paste(ct,"post",sep="_"), 
                                verbose = F, min.pct=0.5,test.use = "wilcox")
      DEGs[[ct]] = DEGs[[ct]][order(abs(DEGs[[ct]]$avg_logFC), decreasing = TRUE), ]
    }, error = function(e) {
      print(e)
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
  return(DEGs)
}

save_DEGs_prior_post <- function(DEGs, summary, f_name) {
  wb <- createWorkbook()
  
  s <- createSheet(wb, "summary")
  addDataFrame(summary, s)

  # DEGs by cell type
  lapply(names(DEGs), function(x) {
    s <- createSheet(wb, x)
    addDataFrame(DEGs[[x]], s)
  })
  
  saveWorkbook(wb, f_name)
}


get_DEGs_summary <- function(DEGs) {
  # filtering only genes over/underexpressed in MeOH samples
  summary = data.frame()
  for (ct in names(DEGs)) {
    # looping through CT
    for (gene in rownames(DEGs[[ct]])) {
      # looping through genes
      if (abs(DEGs[[ct]][gene, "avg_logFC"]) > 0.5) {
        if (gene %in% rownames(summary)) {
          i = if (DEGs[[ct]][gene, "avg_logFC"] > 0.5) 1 else 2
          summary[gene, i] = summary[gene, i] + 1
        } else {
          if (DEGs[[ct]][gene, "avg_logFC"] > 0.5)
            n = data.frame(c(1), c(0))
          else
            n = data.frame(c(0), c(1))
          
          rownames(n) = c(gene)
          summary = rbind(summary, n)
        }
      }
    }
  }
  if (length(names(summary)) > 0) {
    names(summary) = c("greater_than_0.5", "less_than_0.5")
    summary = summary[order(summary[,1] + summary[,2], decreasing = T),] # order by largest to smallest
  } else {
    summary = data.frame("No DEGs > 0.5 || < 0.5" = "None")
  }
  return(summary)
}

