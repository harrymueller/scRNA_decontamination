# source files
source("/data/decont_project/scripts/develop/scripts/general_functions.R")

config = list("methods" = c("soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "decontx:paper", "fastcar", "cellbender"))
load_libraries()

dir = "/data/decont_project/analyses/hgmm12k/results"
output_dir = "/data/decont_project/analyses/hgmm12k/manuscript/scatterplot"

new_names = list("soupx:autoEstCont" = "SoupX with\nautoEstCont",
                 "soupx:background_genes" = "b) SoupX with Three Marker\nGenes per Cell Type",
                 "soupx:top_background_genes" = "SoupX with the\nTop 25 Marker\nGenes",
                 "decontx:no_cell_types" = "DecontX without\nCell Types",
                 "decontx:with_cell_types" = "DecontX with\nCell Types",
                 "decontx:paper" = "d) DecontX with Cell Types\nand `max-iter = 60`",
                 "fastcar" = " \nc) FastCAR",
                 "cellbender" = " \na) CellBender")

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
for (ct in c("mm10", "hg19")) {
  data = data.frame()
  
  # read in data appending to `data`
  for (m in config$methods) {
    temp = read_excel(paste(dir, m, "transcript_origins.xlsx", sep="/"), sheet = ct)
    n = length(rownames(temp))
    to_add = data.frame(rep(new_names[[m]], n), 
                        temp$exo_counts_after,
                        temp$exo_counts_before)
    data = rbind(data, to_add)
  }
  
  names(data) = c("method", "exo_counts_after", "exo_counts_before")
  data["Species"] = ct
  
  dat = rbind(dat, data)
}

# randomise rows of df
set.seed(31415)
rows = sample(nrow(dat))
dat = dat[rows,]
  
for (m in config$methods) {
  all = ggplot(dat[dat$method == new_names[[m]],], aes(x=exo_counts_before, y=exo_counts_after, color=Species)) +
    geom_point(alpha = 0.7) + ggtitle(new_names[[m]]) + 
    geom_abline(slope=1, intercept=0) +
    theme_bw() +
    scale_color_manual(labels = c("Human", "Mouse"), values = c("#5D3A9B", "#E66100"), 
                       guide = guide_legend(override.aes = list(size = 3, alpha = 1) )) + 
    ylab("Post-Decontamination") + xlab("Prior-Decontamination") +
    #ylab("Exogenous UMI Counts Post-Decontamination") + xlab("Exogenous UMI Counts Prior-Decontamination") +
    theme(text=element_text(size=28))

  subset = all + ylim(0, 1000) + xlim(0, 1000)
  all = all + ylim(0, 4000) + xlim(0, 4000)
  
  #subset all + ylim(0,1000) + xlim(0,1000) + 
  #   ggtitle(NULL) +
  #   theme_bw() + theme(axis.line=element_blank(),
  #                           axis.text.x=element_blank(),
  #                           axis.text.y=element_blank(),
  #                           axis.ticks=element_blank(),
  #                           axis.title.x=element_blank(),
  #                           axis.title.y=element_blank(),
  #                           legend.position="none",)
  
  # save plots
  ggsave(paste0(output_dir, "/", m, ".png"), all, width=12, height=10)
  ggsave(paste0(output_dir, "/", m, "_subset.png"), subset, width=12, height=10)
}
