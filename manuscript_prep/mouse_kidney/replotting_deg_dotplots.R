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
OUTPUT_DIR = "/data/decont_project/analyses/mouse_kidney/output_reclus/dotplots"


# create directories
for (dir in c("nine", "specific")) {
  if (!dir.exists(paste(OUTPUT_DIR, dir, sep = "/")))
    dir.create(paste(OUTPUT_DIR, dir, sep = "/"))
}


opts = list(
  s = 14,
  cols = c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000") # use IBM colour-blind palette - https://lospec.com/palette-list/ibm-color-blind-safe
)

# sep = "Under" || "Over"
plot_dotplots_for_n <- function(samples.combined, order_paper, output_dir, DEGs.sep, sep, n, opts, method) {
  Idents(samples.combined) <- "celltype"
  
  # order_paper = CT_ORDER
  # output_dir = OUTPUT_DIR
  # DEGs.sep = DEGs.over
  # sep = "Over"
  # n = 0
  # method = "soupx:autoEstCont"
  
  # choosing genes
  if (n == 9) {
    genes = get_genes_de(DEGs.sep, 8) # > 8
    genes.n = 9
    output_dir = paste0(output_dir, "/nine")
  } else {
    genes = degs_opt_ct(DEGs.sep, 9)
    genes.n = strtoi(genes[1])
    
    if (genes.n == 9) return()
    else print(paste(method, sep, genes.n+1, sep = ", "))
    
    genes = genes[seq(2, length(genes))]
    output_dir = paste0(output_dir, "/specific")
  }
  if (length(genes) == 0) return()
  
  # PLOTTING
  p1 = plot(samples.combined, order_paper, opts, "fresh", genes, "i. Fresh Samples", method)#paste0("DEGs ", sep, "-Expressed in MeOH\n>= ", genes.n, " Celltypes\nFresh Cells"))
  p2 = plot(samples.combined, order_paper, opts, "meoh", genes, "ii. MeOH Samples", method)
  
  #find limits of average expressions
  expr = c(p1$data$avg.exp.scaled, p2$data$avg.exp.scaled)
  expr[is.infinite(expr)] = NA
  lim = c(floor(min(expr, na.rm = T)), ceiling(max(expr, na.rm = T)))
  
  # set legend
  p1 = p1 + scale_color_gradientn(limits = lim, colours=opts$cols) + NoLegend()
  p2 = p2 + scale_color_gradientn(limits = lim, colours=opts$cols) + ylab("")
  p = p1 + p2
  ggsave(paste(output_dir,"/DEG_",method, "_" ,sep,"expr.png",sep=""),
         p,width=16, height=6)
  return(p)
}

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

# loop through number of ct to find num degs between 7 and 10 in more than num_ct cells
degs_opt_ct <- function(DEGs, num_ct) {
  for (i in seq(0, num_ct)) {
    genes <- get_genes_de(DEGs, i)
    if (length(genes) >= 7 && length(genes) <= 11) return(c(i, genes))
  }
  return(c(0, get_genes_de(DEGs, 0))) # default of 1 cell type
}

for (method in config$methods) {
  print(method)
  # method = config$methods[1]
  # load data and ensure metadata is added
  samples.combined = readRDS(paste(INPUT_DIR, method, "Rda", "integrated_rd.Rda", sep = "/"))
  samples.combined = adding_metadata(samples.combined)
  Idents(samples.combined) = "celltype"
  DefaultAssay(samples.combined) <- "RNA"
  
  # # ensure order is correct
  # for (fixation in c("fresh", "meoh")) {
  #   ident = paste0("celltype.", fixation)
  #   samples.combined@meta.data[[ident]] = factor(samples.combined@meta.data[[ident]], levels = c(CT_ORDER, paste0(CT_ORDER, if (fixation == "fresh") "_MeOH" else "fresh")))
  # }
  
  # get precalculated DEGs
  DEGs = readRDS(paste0(INPUT_DIR, "/deg_rda/", method, ".Rda"))
  
  # Seperating all DEGs into overexpressed and underexpressed in MeOH samples
  DEGs.over = get_over_under_DEGs(DEGs, TRUE)   # over-expressed
  DEGs.under = get_over_under_DEGs(DEGs, FALSE) # under-expressed
  
  #DEGs_dotplot_specific_custom(samples.combined, CT_ORDER, OUTPUT_DIR, DEGs.over, DEGs.under)
  plot_dotplots_for_n(samples.combined, CT_ORDER, OUTPUT_DIR, DEGs.over, "Over", 9, opts, method) #MAIN
  plot_dotplots_for_n(samples.combined, CT_ORDER, OUTPUT_DIR, DEGs.under, "Under", 9, opts, method)
  plot_dotplots_for_n(samples.combined, CT_ORDER, OUTPUT_DIR, DEGs.over, "Over", 0, opts, method)
  plot_dotplots_for_n(samples.combined, CT_ORDER, OUTPUT_DIR, DEGs.under, "Under", 0, opts, method)
}
