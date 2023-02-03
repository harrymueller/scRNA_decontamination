# source files
source("/data/decont_project/scripts/develop/scripts/general_functions.R")
source("/data/decont_project/scripts/develop/scripts/mouse_kidney/analysis_functions.R")
source("/data/decont_project/scripts/develop/scripts/mouse_kidney/DEGs.R")

CT_ORDER = rev(c("PT", "Endo", "CD_IC", "MC", "B", "DCT_CNT", "CNT", "DCT", "NK", "T", "CD_PC", "aLOH", "Podo", "MPH", "Fib", "CD_Trans", "Unknown"))


# load libraries
config = list("methods" = c("no_decontamination", "soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "fastcar", "cellbender"), ct_order_dotplots = CT_ORDER)
load_libraries() 

# consts
INPUT_DIR = "/data/decont_project/analyses/mouse_kidney/output_reclus"
OUTPUT_DIR = "/data/decont_project/analyses/mouse_kidney/output_reclus/histograms"

opts = list(
  s = 14,
  cols = c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000") # use IBM colour-blind palette - https://lospec.com/palette-list/ibm-color-blind-safe
)

# create a ggplot obj and return it using the provided parameters
plot <- function(samples.combined, order_paper, opts, fixation, genes, title, method) {
  # combine DCT and CNT into one cell type
  addition = if (fixation == "fresh") "_MeOH" else "_fresh"
  if (method == "cellbender") {
    arr = samples.combined@meta.data[[paste0("celltype.", fixation)]]
    arr[arr == "DCT"] = "DCT_CNT"
    arr[arr == "CNT"] = "DCT_CNT"
    arr[arr == paste0("DCT",addition)] = paste0("DCT_CNT",addition)
    arr[arr == paste0("CNT",addition)] = paste0("DCT_CNT",addition)
    arr[arr == paste0("CD_Trans",addition)] = "Unknown"
    samples.combined@meta.data[[paste0("celltype.", fixation)]] = arr
  }
  Idents(samples.combined) <- paste0("celltype.", fixation)
  order_paper = order_paper[which(order_paper %in% levels(samples.combined))]
  
  # ensure correct order
  samples.combined@active.ident = factor(samples.combined@active.ident, levels = c(order_paper, paste0(order_paper, addition)))
  p = DotPlot(samples.combined, assay = "RNA", idents=order_paper, features = genes, dot.scale = 5, cluster.idents = F, scale.min=0, scale.max=100, scale = TRUE) +
    xlab("DEGs") + ylab("Cell Type") +
    theme(text = element_text(size=opts$s+4), axis.text.y = element_text(size=opts$s), axis.text.x = element_text(size=opts$s,angle = 60, hjust = 1)) +
    ggtitle(title)
  
  return(p)
}

new_names = list("No Decontamination" = "Pre-Decontamination",
                 "SoupX with autoEstCont" = "SoupX with\nautoEstCont",
                 "SoupX with Three Marker Genes per Cell Type" = "SoupX with Three\nMarker Genes per\nCell Type",
                 "SoupX with the Top 25 Marker Genes" = "SoupX with the\nTop 25 Marker\nGenes",
                 "DecontX without Cell Types" = "DecontX without\nCell Types",
                 "DecontX with Cell Types" = "DecontX with\nCell Types",
                 "FastCAR" = "FastCAR",
                 "CellBender" = "CellBender")


degs = read_xlsx(paste(INPUT_DIR, "histograms", "DEGs_Summary.xlsx", sep = "/"))
names(degs) = c("method", "unique", "over", "under")

# method names
for (n in names(new_names)) {
  degs$method[degs$method == n] = new_names[[n]]
}

dat = data.frame()
for (e in c("over", "under")) { #"unique", 
  to_add = degs[c("method", e)]
  names(to_add) = c("method", "expression")
  to_add["type"] = e
  dat = rbind(dat, to_add)
}

dat$method = factor(dat$method, levels = degs$method)

p = ggplot(dat, aes(x = method, y = expression, fill = type)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values = c("#5D3A9B", "#E66100"), labels = c("Higher Expression", "Lower Expression")) +
  xlab("") + ylab("Number of DEGs") + #Decontamination Method
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20)) +
  labs(fill="Differential Expression in\nMethanol-Fixed Samples\nCompared to Fresh Samples")

ggsave(paste(OUTPUT_DIR, "histogram.png", sep="/"), p, width = 15, height = 8)
