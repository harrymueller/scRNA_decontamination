# TODO comment this file

# native_counts = endogenous UMI counts
# non_native = exogenous
# fractions = UMI proportions

identify_transcript_origin <- function (samples.combined) {
  # list of genes and origin
  genes = sapply(rownames(samples.combined), function(x) { str_split(x, "-")[[1]][1] })
  
  # matrix subsets containing only mm || hg genes
  mm_genes = samples.combined@assays$RNA@counts[names(genes)[which(genes == "mm10")],]
  hg_genes = samples.combined@assays$RNA@counts[names(genes)[which(genes == "hg19")],]
  
  # DF (barcode, ct, endo_counts, non_endo_counts, endo_frac)
  transcripts = data.frame("barcode"=names(Idents(samples.combined)), "celltype"=Idents(samples.combined))
  transcripts$barcode = unfactor(transcripts$barcode)
  
  # template
  transcripts$native_counts = rep(NaN, length(transcripts$barcode))
  transcripts$non_native_counts = transcripts$native_counts
  
  # for each subset: for each cell: sum endo transcripts against non-endo transcripts -> calculate fraction of endo
  mm_counts = colSums(as.matrix(mm_genes))
  hg_counts = colSums(as.matrix(hg_genes))
  
  # assigning counts to transcripts DF
  for (i in seq(length(transcripts$barcode))) {
    if (transcripts$celltype[i] == "hg19") {
      transcripts$native_counts[i] = hg_counts[transcripts$barcode[i]]
      transcripts$non_native_counts[i] = mm_counts[transcripts$barcode[i]]
    } else {
      transcripts$native_counts[i] = mm_counts[transcripts$barcode[i]]
      transcripts$non_native_counts[i] = hg_counts[transcripts$barcode[i]]
    }
  }
  
  # native counts over all counts
  transcripts$native_frac = transcripts$native_counts / (transcripts$native_counts + transcripts$non_native_counts)
  transcripts$non_native_frac = 1 - transcripts$native_frac
  
  return(transcripts)
}



# "merge" method transcripts with none transcripts to produce a DF with counts removed etc.
combine_transcripts <- function (transcripts_none, transcripts_method) {
  transcripts = data.frame(barcodes = transcripts_none$barcode, celltype = transcripts_none$celltype)
  transcripts$barcodes = unfactor(transcripts$barcodes)
  transcripts_method = transcripts_method[order(match(transcripts_method$barcode, transcripts_none$barcode)),]


  loop = c("endo_counts" = "native_counts",
           "exo_counts"  = "non_native_counts",
           "endo_fract"  = "native_frac",
           "exo_fract"   = "non_native_frac")

  for (n in names(loop)) {
    transcripts[paste(n, "_before", sep="")] = transcripts_none[[loop[[n]]]]
    transcripts[paste(n, "_after", sep="")] = transcripts_method[[loop[[n]]]]
  }

  # count differences:
  transcripts["endo_counts_diff"] = transcripts["endo_counts_after"] - transcripts["endo_counts_before"]
  transcripts["exo_counts_diff"] = transcripts["exo_counts_after"] - transcripts["exo_counts_before"]
  
  # fraction differences
  # OLD method: (frac_after - frac_before) / frac_before
  transcripts["endo_fract_diff"] = transcripts["endo_counts_diff"] / transcripts["endo_counts_before"]
  transcripts["exo_fract_diff"] = transcripts["exo_counts_diff"] / transcripts["exo_counts_before"]

  return(transcripts)
}



summarise_transcripts <- function (transcripts_combined) {
  # template DF
  summ = data.frame("celltypes" = rep(c("hg19", "mm10", "TOTAL"), 3), 
                    "measure"=c(rep("mean", 3), rep("median",3), rep("std", 3)),
                    "exo_fract_after"=rep(NaN,9), "exo_fract_diff"=rep(NaN,9), 
                    "endo_counts_diff"=rep(NaN,9), "endo_fract_diff"=rep(NaN,9))
  
  summ$celltypes = unfactor(summ$celltypes)
  summ$measure = unfactor(summ$measure)
  
  # loop through each "ct" (inc. overall)
  for (ct in summ$celltypes[1:3]) {
    # create subsets of transcripts
    if (ct != "TOTAL") {
      transcripts_subset = transcripts_combined[transcripts_combined$celltype == ct,]
      i = if (ct == "hg19") 1 else 2 # index for summary_transcripts
    } else {
      transcripts_subset = transcripts_combined
      i = 3
    }
    
    # mean, median and standard deviation for each
    for (n in c("exo_fract_after", "exo_fract_diff", "endo_counts_diff", "endo_fract_diff")) {
      summ[i, n]   = mean(transcripts_subset[[n]])
      summ[i+3, n] = median(transcripts_subset[[n]])
      summ[i+6, n] = sd(transcripts_subset[[n]])
    }
  }
  
  return(summ)
}



# summarise and save the measure of dataset prior to decontamination
# similar to summarise_transripts & save_summary_transcripts
summarise_transcripts_before_decont = function (transcripts) {
  # TODO add median
  # template DF
  summ = data.frame("celltypes" = rep(c("hg19", "mm10", "TOTAL"), 3), 
                    "measure"=c(rep("mean", 3), rep("median", 3), rep("std",3)),
                    "native_counts"=rep(NaN,9), "non_native_counts"=rep(NaN,9),                       
                    "native_frac"=rep(NaN,9), "non_native_frac"=rep(NaN,9))
  
  summ$celltypes = unfactor(summ$celltypes)
  summ$measure = unfactor(summ$measure)

  for (ct in summ$celltypes[1:3]) {
    # create subsets of transcripts
    if (ct != "TOTAL") {
      transcripts_subset = transcripts[transcripts$celltype == ct,]
      i = if (ct == "hg19") 1 else 2 # index for summary_transcripts
    } else {
      transcripts_subset = transcripts
      i = 3
    }

    # mean, median and standard deviation for each
    for (n in c("native_counts", "non_native_counts", "native_frac", "non_native_frac")) {
      summ[i, n]   = mean(transcripts_subset[[n]])
      summ[i+3, n] = median(transcripts_subset[[n]])
      summ[i+6, n] = sd(transcripts_subset[[n]])
    }
  }

  # round and convert fractions to percentages
  summ[3:4] = lapply(summ[3:4], round, 2)
  summ[5:6] = lapply(summ[5:6], function(x) {
    round(x, 5)*100
  })

  new = data.frame("Contamination Fraction Statistics"=rep("", 4),
                   "celltype" = c("measure", "hg19", "mm10", "all"))
  

  # reformat summ for better output
  loop = c("endogenous count" = "native_counts",
           "exogenous count" = "non_native_counts",
           "endogenous fraction" = "native_frac",
           "exogenous fraction" = "non_native_frac")

  for (n in names(loop)) {
    new[n]            = c("mean", abs(summ[1:3, loop[[n]]]))
    new[paste(n,"")]  = c("median", abs(summ[4:6, loop[[n]]]))
    new[paste(n," ")] = c("std", abs(summ[7:9, loop[[n]]]))
  }
  
  # save to excel: first sheet is summary, then full DF (separated by CT)
  wb <- createWorkbook()
    
  s <- createSheet(wb, "summary")
  addDataFrame(t(new), s, row.names = TRUE, col.names = FALSE)
  
  for (ct in c("hg19", "mm10")) {
    s = createSheet(wb, ct)
    addDataFrame(transcripts[transcripts$celltype == ct,], s, row.names = FALSE, col.names = TRUE)
  }
  
  saveWorkbook(wb, paste(files$output, "/transcript_origins.xlsx",sep="/"))
}



# Given the summary and transcripts DF -> saves to an xlsx file 
save_summary_transcripts <- function (transcripts, summ) {
  # round to 5 sig fig and convert to percentages
  summ[c(3,4,6)] = lapply(summ[c(3,4,6)], function(x) {
    round(x, 5)*100
  })
  summ[5] = lapply(summ[5], round, 2)

  # template DF for output
  new = data.frame("Contamination Fraction Statistics"=rep("", 4),
                   "celltype" = c("measure", "hg19", "mm10", "all"))
  
  # reformat summ for better output
  loop = c("exogenous percentage removed" = "exo_fract_diff",
           "exogenous percentage remaining" = "exo_fract_after",
           "endogenous counts removed" = "endo_counts_diff",
           "endogenous percentage removed" = "endo_fract_diff")

  for (n in names(loop)) {
    new[n]            = c("mean", abs(summ[1:3, loop[[n]]]))
    new[paste(n,"")]  = c("median", abs(summ[4:6, loop[[n]]]))
    new[paste(n," ")] = c("std", abs(summ[7:9, loop[[n]]]))
  }

  # save to excel: first sheet is summary, then full DF (separated by CT)
  wb <- createWorkbook()
    
  s <- createSheet(wb, "summary")
  addDataFrame(t(new), s, row.names = TRUE, col.names = FALSE)
  
  # adding "raw" data
  for (ct in c("hg19", "mm10")) {
    s = createSheet(wb, ct)
    addDataFrame(transcripts[transcripts$celltype == ct,], s, row.names = FALSE, col.names = TRUE)
  }
  
  saveWorkbook(wb, paste(files$output, "/transcript_origins.xlsx",sep="/"))
}



# creates a plot of all barcodes, before and after decont, against mouse transcripts and human transcripts
plot_transcripts <- function (transcripts) {
  for (i in c("before", "after")) {
    # select certain columns from transcripts and rename
    subset = transcripts[c("barcodes", "celltype", paste("endo_counts", i, sep="_"), paste("exo_counts", i, sep="_"))]
    names(subset) = c("barcodes", "celltype", "human_counts", "mouse_counts")
    
    # swap exogenous and endogenous UMIs for mouse cells - so transcripts are specific to either human or mouse 
    # instead of generic endogenous vs exogenous
    temp = subset$human_counts[subset$celltype=="mm10"]
    subset$human_counts[subset$celltype=="mm10"] = subset$mouse_counts[subset$celltype=="mm10"]
    subset$mouse_counts[subset$celltype=="mm10"] = temp
    
    # add a group column (for plot labels)
    subset["group"] = sapply(subset$celltype, function(ct) {
      paste(ct, if (i == "before") "prior-decont" else "post-decont", sep=" ")
    },USE.NAMES=F)
    
    if (i == "before")
      plot_df = subset
    else
      plot_df = rbind(plot_df, subset)
  }

  # plot of both cell types (x-axis = human UMIs & y-axis = mouse UMIs)
  all <- ggplot(plot_df, aes(x=human_counts, y=mouse_counts, color=group)) + 
                  geom_point() + ggtitle("UMI Counts of All Cells") +
                  ylab("Mouse UMI Counts") + xlab("Human UMI Counts") +
                  scale_color_manual(values=c("#27AE60", "#8E44AD", "#E67E22","#3498DB"), name="Species and Prior- or\nPost-Decontamination") +
                  theme(text=element_text(size=16, family="TT Times New Roman"))

  # mouse UMI counts in any human cells with mouse counts > 1000 are set to 1000 
  plot_df[plot_df$celltype == "hg19" & plot_df$mouse_counts > 1000,"mouse_counts"] = 1000

  # plot of ONLY human cells (with mouse UMIs <= 1000) - will remove some cells
  human = ggplot(plot_df[plot_df$celltype == "hg19",], aes(x=human_counts, y=mouse_counts, color=group)) +
                  geom_point() + ggtitle("UMI Counts of Human Cells") + 
                  ylab("Mouse UMI Counts (0-1000)") + xlab("Human UMI Counts") + ylim(0, 1000) +
                  scale_color_manual(values=c("#27AE60", "#8E44AD"), name="Species and Prior- or\nPost-Decontamination") + 
                  theme(text=element_text(size=16, family="TT Times New Roman"))

  # human UMI counts in any mouse cells with human counts > 1000 are set to 1000
  plot_df[plot_df$celltype == "mm10" & plot_df$human_counts > 1000,"human_counts"] = 1000

  # plot of ONLY mouse cells (with human UMIs <= 1000) - will remove some cells
  mouse = ggplot(plot_df[plot_df$celltype == "mm10",], aes(x=mouse_counts, y=human_counts, color=group)) + 
                  geom_point() + ggtitle("UMI Counts of Mouse Cells") + 
                  ylab("Human UMI Counts (0-1000)") + xlab("Mouse UMI Counts") + ylim(0, 1000) +
                  scale_color_manual(values=c("#E67E22", "#3498DB"), name="Species and Prior- or\nPost-Decontamination") +
                  theme(text=element_text(size=16, family="TT Times New Roman"))

  # save plots
  ggsave(paste(files$output, "plots/UMI_plot_for_all_cells.png", sep="/"), all, width=12, height=10)
  ggsave(paste(files$output, "plots/UMI_plots_for_each_celltype.png", sep="/"), (human / mouse), width=12, height=10)
  ggsave(paste(files$output, "plots/UMI_plots_combined.png", sep="/"), all / (human / mouse), width=12, height=20)
}