# source files
source("/data/decont_project/scripts/develop/scripts/general_functions.R")

config = list("methods" = c("soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "decontx:paper", "fastcar", "cellbender"))
load_libraries()

dir = "/data/decont_project/analyses/hgmm12k/results"
output_dir = "/data/decont_project/analyses/hgmm12k/manuscript/violin"

new_names = list("soupx:autoEstCont" = "SoupX with\nautoEstCont",
                 "soupx:background_genes" = "SoupX with Three\nMarker Genes per\nCell Type",
                 "soupx:top_background_genes" = "SoupX with the\nTop 25 Marker\nGenes",
                 "decontx:no_cell_types" = "DecontX without\nCell Types",
                 "decontx:with_cell_types" = "DecontX with\nCell Types",
                 "decontx:paper" = "DecontX with\nCell Types and\n`max-iter = 60`",
                 "fastcar" = "FastCAR",
                 "cellbender" = "CellBender")

# new_names = list("soupx:autoEstCont" = "SoupX 1",
#                  "soupx:background_genes" = "SoupX 2",
#                  "soupx:top_background_genes" = "SoupX 3",
#                  "decontx:no_cell_types" = "DecontX 1",
#                  "decontx:with_cell_types" = "DecontX 2",
#                  "decontx:paper" = "DecontX 3",
#                  "fastcar" = "FastCAR",
#                  "cellbender" = "CellBender")

dat = NULL
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
  data["Species"] = ct

  dat = rbind(dat, data)
}

# PLOTTING
col = "exo_removed"#for (col in c("exo_removed")){#, "exo_remaining", "endo_removed")) {
dat["Y"] = dat[col]

p = ggplot(dat, aes(x = method, y = Y, fill = Species)) + 
  geom_violin(width = 1, color = alpha('black', 0.5)) +
  scale_fill_manual(values = c("#5D3A9B", "#E66100"), labels = c("Human", "Mouse")) +
  xlab("") + ylab("Exogenous UMIs Removed (%)") + #Decontamination Method
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20))
    
ggsave(paste(output_dir, "/violin_plot_", col, ".pdf", sep=""), p, width = 12, height = 6.5)
ggsave(paste(output_dir, "/violin_plot_", col, ".png", sep=""), p, width = 12, height = 7)
