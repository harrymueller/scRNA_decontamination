# source files
source("scripts/general_functions.R")

config = list("methods" = c("soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "fastcar", "cellbender"))
load_libraries()

dir = "/data/Perkins/hgmm12k/results"
output_dir = "/data/Perkins/hgmm12k/results/violin_plots"

# read in data
for (ct in c("hg19", "mm10")) {
  data = data.frame()
  
  # read in data appending to `data`
  for (m in config$methods) {
    temp = read_excel(paste(dir, m, "transcript_origins.xlsx", sep="/"), sheet = ct)
    n = length(rownames(temp))
    to_add = data.frame(rep(m, n), 
                        (temp$exo_counts_before - temp$exo_counts_after) / temp$exo_counts_before * 100,
                        temp$exo_counts_after / (temp$exo_counts_after + temp$endo_counts_after) * 100,
                        (temp$endo_counts_before - temp$endo_counts_after) / temp$endo_counts_before * 100)
    data = rbind(data, to_add)
  }
  
  names(data) = c("method", "exo_removed", "exo_remaining", "endo_removed")
  
  # plot
  for (col in c("exo_removed", "exo_remaining", "endo_removed")) {
    data["Y"] = data[col]
    p = ggplot(data, aes(x = method, y = Y)) + 
      geom_violin(trim = T) +
      ggtitle(paste("Violin Plots of", col)) +
      ylab(paste(col, "%")) + #scale_y_continuous(limits=c(0, 3)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text=element_text(size=16, family="TT Times New Roman"))
    ggsave(paste(output_dir, "/violin_plot_", ct, "_", m, "_", col, ".png", sep=""), p, width = 9, height = 9)
  }
}

# format data for plotting


# plot