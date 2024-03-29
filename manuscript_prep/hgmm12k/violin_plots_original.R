# source files
source("/data/decont_project/scripts/develop/scripts/general_functions.R")

config = list("methods" = c("soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "decontx:paper", "fastcar", "cellbender"))
load_libraries()

dir = "/data/decont_project/analyses/hgmm12k/results"
output_dir = "/data/decont_project/analyses/hgmm12k/manuscript/violin"

new_names = list("soupx:autoEstCont" = "SoupX w/\nAutoEstCont",
				 "soupx:background_genes" = "SoupX w/ 3 Marker\nGenes per CT",
				 "soupx:top_background_genes" = "SoupX w/ Top 25\nMarkers",
				 "decontx:no_cell_types" = "DecontX without\nCell Types",
				 "decontx:with_cell_types" = "DecontX with\nCell Types",
				 "decontx:paper" = "DecontX with\nParameters from\nthe Publication",
				 "fastcar" = "FastCAR",
				 "cellbender" = "CellBender")

new_names = list("soupx:autoEstCont" = "SoupX 1",
                 "soupx:background_genes" = "SoupX 2",
                 "soupx:top_background_genes" = "SoupX 3",
                 "decontx:no_cell_types" = "DecontX 1",
                 "decontx:with_cell_types" = "DecontX 2",
                 "decontx:paper" = "DecontX 3",
                 "fastcar" = "FastCAR",
                 "cellbender" = "CellBender")

# read in data
for (ct in c("hg19", "mm10")) {
  data = data.frame()
  
  # read in data appending to `data`
  for (m in config$methods) {
    temp = read_excel(paste(dir, m, "transcript_origins.xlsx", sep="/"), sheet = ct)
    n = length(rownames(temp))
    to_add = data.frame(rep(new_names[[m]], n), 
                        (temp$exo_counts_before - temp$exo_counts_after) / temp$exo_counts_before * 100,
                        temp$exo_counts_after / (temp$exo_counts_after + temp$endo_counts_after) * 100,
                        (temp$endo_counts_before - temp$endo_counts_after) / temp$endo_counts_before * 100)
    data = rbind(data, to_add)
  }
  
  names(data) = c("method", "exo_removed", "exo_remaining", "endo_removed")
  
  # plot
  for (col in c("exo_removed")){#, "exo_remaining", "endo_removed")) {
    data["Y"] = data[col]
    p = ggplot(data, aes(x = method, y = Y)) + 
      geom_violin(width = 1, fill = alpha(if (ct == "hg19") '#5D3A9B' else "#E66100", 1.0), color = alpha('black', 0.5)) +
      #ggtitle(paste("Proportion of Exogenous UMIs Removed from", if (ct == "hg19") "Human" else "Mouse" ,"\nCells for each Decontamination Method")) +
      ggtitle(paste("", if (ct == "hg19") "A. Human" else "B. Mouse" ,"Cells")) +
  	  xlab("") + 
  	  ylab("Exogenous UMIs Removed (%)") + 
  	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20))
    ggsave(paste(output_dir, "/violin_plot_", ct, "_", col, ".png", sep=""), p, width = 7, height = 7)
  }
}
