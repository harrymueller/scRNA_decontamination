# Functions relating to the analysis / plotting



################################################################################################
### Adds some metadata to the combined seurat object
################################################################################################
adding_metadata <- function(samples.combined) {
  # assumes Idents is celltype and contains $preservation
  if (!("celltype" %in% colnames(samples.combined@meta.data))) {
    samples.combined@meta.data$celltype = Idents(samples.combined)
    # celltype_method
    samples.combined@meta.data$celltype_method = paste(samples.combined$celltype, samples.combined$preservation, sep="_")
  }
  
  if (F) {
    # LEGACY CODE
    samples.combined$celltype_method <- Idents(samples.combined)
    
    samples.combined$celltype <- sapply(samples.combined$celltype_method, function(annotation) {
      split = c(str_split(annotation, "_", simplify = T))
      return(paste(split[1:length(split)-1], collapse="_"))
    }, simplify = TRUE)
  }
  
  Idents(samples.combined) = "celltype"
  
  # changing ct annotations of CD_Trans to unknown
  if (any(Idents(samples.combined)=="CD_Trans"))
    Idents(samples.combined)[which(Idents(samples.combined)=="CD_Trans")] = "Unknown"
  
  order_paper = config$ct_order_dotplots[which(config$ct_order_dotplots %in% levels(samples.combined))]
  
  #samples.combined <- readRDS(paste(output, "samples_integrated_rd.Rda", sep="/"))
  DefaultAssay(samples.combined) <- "RNA"
  
  if (!("celltype.fresh" %in% colnames(samples.combined@meta.data))) {
    Idents(samples.combined) = samples.combined$celltype_method
    # Following relates to plotting the dotplots
    samples.combined$celltype.fresh = unlist(lapply(samples.combined$celltype_method, function(x) {
      end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="fresh")
      x <- if (end) str_sub(x,end=-7) else x
    }))
    samples.combined$celltype.meoh = unlist(lapply(samples.combined$celltype_method, function(x) {
      end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="MeOH")
      x <- if (end) str_sub(x,end=-6) else x
    }))}
  return(samples.combined)
}


################################################################################################
# Plots of differentially expressed genes
################################################################################################
analyse_DEGs <- function(samples.combined) {
  Idents(samples.combined) <- "celltype_method"

  # Getting DEGs
  DEGs = get_DEGs(samples.combined, paste(files$output,"/Rda/DEGs_per_celltype.Rda",sep=""))

  # Seperating all DEGs into overexpressed and underexpressed
  DEGs.over = get_over_under_DEGs(DEGs, TRUE)   #DEGs, bool_whether_overexpressed
  DEGs.under = get_over_under_DEGs(DEGs, FALSE) #DEGs, bool_whether_overexpressed

  # Identify which DEGs are over/under expressed in at least 9 cell types
  genes.over <- get_genes_de(DEGs.over, 8)   #DEGs.selected, num_cells
  genes.under <- get_genes_de(DEGs.under, 8) #DEGs.selected, num_cells

  ### Saving DEGs to excel file
  #save_DEGs(DEGs, paste(files$output,"/DEGs.xlsx",sep=""), genes.over, genes.under) #DEGs, f_name, genes.over, genes.under
  
  ## Plotting DEGs
  print("Plotting DEGs")
  #DEGs_histogram(DEGs) #DEGs, f_name
  DEGs_dotplot_specific(samples.combined,config$ct_order_dotplots, paste(files$output, "plots", sep="/"), DEGs.over, DEGs.under)
  #DEGs_dotplot_over_under_expression(samples.combined, paste(files$output,"/plots/DEG_higher_expression_9_celltypes.png",sep=""), 
  #                                   config$ct_order_dotplots, config$genes_ct_dotplots,
  #                                   genes.over, genes.under)
}



################################################################################################
# Plots related to reclustering
################################################################################################
analyse_recluster <- function(samples.combined, current_method) {
  print("Plotting CT after reclustering")
  # save & get contingency table of ct numbers
  tables = contigency_table_ct(config$sample_ids, files$CellAnnotations, 
                               paste(files$output, "/new_clus.tsv",sep=""), 
                               paste(files$output, "/clus_contingency_table.xlsx",sep=""),
                               config$is_xlsx[[current_method]])
  
  # plotting module score differences
  plot_module_score_hist(samples.combined, tables, paste(files$output, "/plots/module_score_diff.png",sep=""))

  # barplot of ct changes
  barplots_ct_props(tables$fresh$table, paste(files$output, "/plots/barplot_ct_prop_fresh.png",sep=""))
  barplots_ct_props(tables$meoh$table, paste(files$output, "/plots/barplot_ct_prop_meoh.png",sep=""))
}






################################################################################################
# Plots a pie chart showing the relative proportions of cell-types within all samples
# Plots a pie chart for each preservation method
################################################################################################
plot_pie_ct <- function (samples, method, ident_name) {
  output = files$OcraRelDir
  plot_cts = config$pie_plot_cts

  # plotting pie charts for each preservation technique
  Idents(samples) = ident_name
  samples <- SplitObject(samples)
  plot_cts = plot_cts[which(plot_cts %in% unique(samples[[1]]@meta.data$celltype))]
	
  lapply(samples, function(x) {
    type=unique(x@meta.data[[ident_name]])
    if (type != "decont") {
      df = data.frame(xtabs(~x@meta.data$celltype))
      
      df$labels = unfactor(df[,1])
      df[,1]=NULL

      # creates a new category for other, containing all cts not in plot_cts
      df = rbind(df[which(df$labels %in% plot_cts),], c(sum(df$Freq[which(!(df$labels %in% plot_cts))]),"Other"))
      
      df$Freq = as.integer(df$Freq)
      
      # proportions
      df$prop = df$Freq / sum(df$Freq) *100
      
      # ordering cts - for comparing different analyses
      df = df[order(match(df$labels, c(plot_cts, "Other"))),]

      df$labels = c(plot_cts, "Other")
      
      # Labels 
      df$labels = factor(df$labels, levels=df$labels)

      # plotting
      fig <- plot_ly(df, labels = ~labels, values = ~prop, type = 'pie',textinfo = 'label+percent', sort=F, textfont = list(size = 20))
      fig <- fig %>% layout(#title = paste(method, " (", type,"); Cell type proportions", sep=""),
                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            margin=list(l=50, r=50, t=150, b=50),
                            showlegend = F)

      # using ocra to save as image
      # always displays error - even though it saves
      library("processx")
      tryCatch(expr = {
        # orca would only work on my PC by turning on debug mode - then closing the HTML popup
        orca(fig, file=paste(output, "/plots/", ident_name, "_", method, "_", type, "_ct_pie.png",sep=""),scale=3, verbose = TRUE, debug = T)
      })
    }
  })
}



################################################################################################
# Saving some information on the recluster to a xlsx document
# inc. contigency tables (fresh & meoh), ari, nmi, and difference of unknowns (orig to this analysis)* not strictly accurate
################################################################################################
contigency_table_ct <- function (all_sample_ids, cell_annotations_path, new_clus_path, output_file, annotations_is_xlsx) {
  ### Getting contigency tables (and cts)
  fresh = get_cont_table(all_sample_ids[1:3], cell_annotations_path, new_clus_path, annotations_is_xlsx)
  meoh = get_cont_table(all_sample_ids[4:6], cell_annotations_path, new_clus_path, annotations_is_xlsx)
  
  ### Dataframe containing ARI / NMI / Change in num unknowns -> for adding to workbook
  stat_df = data.frame(
    "Data"= c("Fresh", "MeOH"),
    "ARI" = c(aricode::ARI(fresh$o_ct, fresh$r_ct),aricode::ARI(meoh$o_ct, meoh$r_ct)),
    "NMI" = c(aricode::NMI(fresh$o_ct, fresh$r_ct),aricode::NMI(meoh$o_ct, meoh$r_ct)),
    "Change in Unknown Cells" = c(
      sum(fresh$r_ct == "Unknown") - sum(fresh$o_ct == "Unknown"),
      sum(meoh$r_ct == "Unknown") - sum(meoh$o_ct == "Unknown")
    )
  )
  
  wb <- createWorkbook()
  # adding ari and nmi
  s <- createSheet(wb, "ari_nmi")
  addDataFrame(t(stat_df), s, row.names = T, col.names = F)
  
  # adding fresh contigency table
  notes_df = data.frame(
    "1"="*Cell types on the top are the cell types after reclustering post-decontamination",
    "2"="*Cell types on the left are the original cell types",
    "3"="Changed% Column is the percentage of cells that changed annotation from the given cell type relative to the starting number of cells within that cell type.",
    "4"="Changed% Row is the percentage of cells that changed annotation to the given cell type relative to the number of cells that kept the same annotation for the cell type.",
    "5"="Changed% x Changed% is the percentage of cells that changed annotation relative to the total count of cells"
  )

  s <- createSheet(wb, "Fresh")
  addDataFrame(
    t(notes_df),s, row.names=F, col.names=F
  )

  f = add_percent_changed_clus_cont(fresh)
  addDataFrame(f, s, row.names = T, col.names = T,startRow = 7)
  
  # adding meoh contigency table
  s <- createSheet(wb, "MeOH")
  addDataFrame(
    t(notes_df),s, row.names=F, col.names=F
  )
  f = add_percent_changed_clus_cont(meoh)
  addDataFrame(f, s, row.names = T, col.names = T,startRow = 7)
  
  # save workbook
  saveWorkbook(wb, output_file)
  
  return(list(fresh=fresh, meoh=meoh))
}

add_percent_changed_clus_cont <- function (tab) {
  data = data.frame(rbind(tab$tab))
  l = length(names(data))-1

  data[["Changed%"]] = rep(0,l+1)
  data = rbind(data, rep(0, l+1))
  rownames(data)[l+2] = "Changed%"

  total_changed = 0

  for (ct in names(data)[seq(l)]) {
    ct_new = data["Sum", rownames(data) == ct]
    ct_stayed = data[ct, rownames(data) == ct]
    ct_orig = data$Sum[rownames(data) == ct]

    data[["Changed%"]][rownames(data) == ct] = paste0(round((ct_orig - ct_stayed) / ct_orig,3),"%")
    data["Changed%", rownames(data) == ct] = paste0(round((ct_new - ct_stayed) / ct_stayed,3),"%")

    total_changed = total_changed + (ct_orig - ct_stayed)
  }

  data["Changed%","Changed%"] = paste0(round(total_changed / data["Sum", "Sum"],3),"%")

  return(data)
}


################################################################################################
# Returns a list containing a contigency table, original cell types, and new cell types
################################################################################################
get_cont_table <- function (all_sample_ids, cell_annotations_path, new_clus_path, annotations_is_xlsx) {
  orig_cts = data.frame()
  # getting all barcodes and cts
  for (id in all_sample_ids) {
    t = get_clusters(cell_annotations_path, id, annotations_is_xlsx)
    orig_cts = rbind(orig_cts, data.frame("barcodes"=names(t), "ct"=t, "sample_id"=id))
  }
  
  # formatting data frame
  orig_cts$barcodes = unfactor(orig_cts$barcodes)
  orig_cts$id = paste(unfactor(orig_cts$sample_id), orig_cts$barcodes, sep="_")
  
  orig_cts$ct = unfactor(orig_cts$ct)
  orig_cts$ct[which(orig_cts$ct=="DCT")]="DCT_CNT"
  orig_cts$ct[which(orig_cts$ct=="CNT")]="DCT_CNT"
  orig_cts$ct = factor(orig_cts$ct,levels=unique(orig_cts$ct))
  
  orig_cts[,c("barcodes", "sample_id")]=NULL
  
  ## Getting clusters after reclustering
  new_cts = read.table(new_clus_path,header = F, sep="\t",skip = 1)
  names(new_cts) = c("id","ct")
  new_cts$id = unfactor(new_cts$id)
  
  # Ordering new table by cell barcodes of original
  new_cts = new_cts[which(new_cts$id %in% orig_cts$id),]
  orig_cts = orig_cts[which(orig_cts$id %in% new_cts$id),]
  new_cts=new_cts[order(match(new_cts$id, orig_cts$id)),]
  
  OriginalCellTypes = orig_cts$ct
  names(OriginalCellTypes) = orig_cts$id
  ReclusteredCellTypes = new_cts$ct
  
  ReclusteredCellTypes = factor(unfactor(ReclusteredCellTypes), levels=levels(OriginalCellTypes))
  names(ReclusteredCellTypes) = new_cts$id
  tab = addmargins(xtabs(~ OriginalCellTypes + ReclusteredCellTypes,))

  return(list(
    table=tab, o_ct=OriginalCellTypes, r_ct=ReclusteredCellTypes
  ))
}

################################################################################################
# plot a density histogram of the differences in module score between two highest scores
################################################################################################
plot_module_score_hist <- function (seurat, tables, output) {
  plots = list()
  for (m in c(1,2)) {
    t = tables[[m]]
    changed_barcodes = names(t$o_ct)[which(t$o_ct != t$r_ct)]
    
    values = seurat@meta.data$module_score_diff[which(colnames(seurat)%in% changed_barcodes)]
    
    df = data.frame("barcodes"=changed_barcodes, "values"=values)
    
    titles = c("Fresh", "MeOH")
    plots[[m]] = ggplot(aes(x=values),data = df) + stat_density() + ggtitle(paste(titles[[m]]," (",length(values)," cells changed)", sep="")) +
      xlab(if (m == 1) NULL else "Difference in Module Scores") + ylab("Density")

    if (config$fonts) plots[[m]] = plots[[m]] + theme(text=element_text(size=16, family="TT Times New Roman"))
  }
  
  # getting maximum values (to set axis limits to be equal)
  x_max = max(ggplot_build(plots[[1]])$layout$panel_scales_x[[1]]$range$range, 
              ggplot_build(plots[[2]])$layout$panel_scales_x[[1]]$range$range)
  y_max = ceiling(max(ggplot_build(plots[[1]])$layout$panel_scales_y[[1]]$range$range, 
                      ggplot_build(plots[[2]])$layout$panel_scales_y[[1]]$range$range))
  
  plots[[3]] = plots[[1]] / plots[[2]] + 
    plot_annotation(title="Density plots of differences in module score for cells that changed annotation") & 
    xlim(0, x_max) & ylim(0,y_max)
  
  ggsave(output,plots[[3]],width=7,height=6)
}

################################################################################################
# Given a contigency table, plots the change in ct proportion (original and new) for each cell type
################################################################################################
barplots_ct_props <- function(table, output) {
  t = data.frame(table)

  # only keep cts that both have
  n_cts = unique(t$ReclusteredCellTypes)
  t = t[which(t$OriginalCellTypes %in% n_cts),]	
  
  # (length-1) removes sum entry
  n = length(n_cts)-1
  orig_ct_sums = t$Freq[which(t$ReclusteredCellTypes=="Sum")][1:n]
  new_ct_sums = t$Freq[which(t$OriginalCellTypes=="Sum")][1:n]
  
  # removing sum columns
  t = t[which(t$OriginalCellTypes!="Sum"),]
  t = t[which(t$ReclusteredCellTypes!="Sum"),]
  
  # number of cells that remained the cell type
  same_ct = t$Freq[which(t$OriginalCellTypes == t$ReclusteredCellTypes)]
  
  # number of cells changing to a ct and from a ct
  changed_to_ct = new_ct_sums-same_ct
  changed_from_ct = orig_ct_sums - same_ct
  
  # 
  #cbind(levels(t$OriginalCellTypes)[1:n],orig_ct_sums, new_ct_sums, same_ct, changed_to_ct, changed_from_ct)
  
  # combining into 1 df
  props = data.frame("CT_Same_Orig"=same_ct/orig_ct_sums, "CT_changed_from_Original"=changed_from_ct/orig_ct_sums,
                     "CT_Same_New"=same_ct/new_ct_sums, "CT_changed_to_New"=changed_to_ct/new_ct_sums,    
                     "labels"=unique(t$OriginalCellTypes)[1:n])
  
  # "melting" the df, for plotting
  melted = melt(props, "labels")
  
  melted$CellTypes = ""
  melted[melted$variable == "CT_Same_Orig",]$CellTypes<-"1. Original CT"
  melted[melted$variable == "CT_changed_from_Original",]$CellTypes<-"1. Original CT"
  
  melted[melted$variable == "CT_Same_New",]$CellTypes<-"2. New CT"
  melted[melted$variable == "CT_changed_to_New",]$CellTypes<-"2. New CT"
  melted$variable = factor(unfactor(melted$variable), levels=rev(c("CT_Same_New", "CT_changed_to_New","CT_Same_Orig", "CT_changed_from_Original")))
  
  melted$num = ""
  
  l = length(table[,1])
  
  to_add = table[1:n,l]
  melted[melted$variable=="CT_changed_from_Original",]$num = to_add[order(match(names(to_add),melted[melted$variable=="CT_changed_from_Original",]$labels))]
  to_add = table[l,1:n]
  melted[melted$variable=="CT_changed_to_New",]$num = to_add[order(match(names(to_add),melted[melted$variable=="CT_changed_to_New",]$labels))]
  
  # replacing NaN w/ 0
  melted$value[which(is.nan(melted$value))]=0

  # plotting
  p <- ggplot(melted, aes(x = CellTypes, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    facet_grid(~ labels) + 
    scale_fill_manual(values = c("red3", "blue3", "orangered1","royalblue2"), 
                      name = "Proportion of barcodes given a\nCT (relative to bar label)",
                      labels=c("Changed to a different CT", "Stayed the same CT","Changed from a different CT", "Stayed the same CT")) +
    geom_text(aes(label=num), size = 3, hjust = 0.5, vjust = 2, position = "stack",color="white") +
    xlab("Cell Types (and proportions origin)") + ylab("Proportion of Cell Barcodes given a CT")

  if (config$fonts)
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, family="TT Times New Roman")) 
  else
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # saving plot
  ggsave(output,p,width=14,height=7)
}