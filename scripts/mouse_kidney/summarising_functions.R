# Summary functions




################################################################################################
# Loads all required libraries
################################################################################################
deg_summary <- function () {
  # variables
  input_path = config$output_dir
  output_path = c(paste(config$output_dir, "summary", "Summary_Histogram.png", sep="/"),
                  paste(config$output_dir, "summary", "DEGs_Summary.xlsx", sep="/"))
  methods = config$methods
  labels = config$summary_histogram_labels
  n = length(methods)
  
  # creating files df
  files = data.frame("filenames" = rep(1, n))
  files$filenames = sapply(methods, function(x) {
    return(paste(input_path, x, "DEGs.xlsx", sep="/"))
  }, USE.NAMES = FALSE)
  
  files$labels = labels
  
  # creating df for counts
  counts = data.frame(nums = rep(0,2*n), method=c(files$labels,files$labels), expression=c(rep("Over-Expressed",n),rep("Under-Expressed",n)))

  unique_degs = list("counts" = rep(-1,n))
  
  # looping through files and sheets to get counts of DEGs and unique DEGs
  suppressMessages({
    for (i in seq(length(files$filenames))) {
      file = files$filenames[i]
      label = files$labels[i]
      
      #iles_index = which(files$filenames==file)
      all_degs = c()

      sheets = excel_sheets(file)
      range = seq(1, length(sheets))

      range = range[which(!(sheets %in% c("overexpressed_genes_>8_celltype", "CD_Trans", "Unknown")))] 
      
      # looping through valid sheets
      for (j in range) {
        DEGs = read_excel(file, sheet=j)

        if (length(DEGs) > 0) {
          all_degs = append(all_degs, DEGs[[1]])
          counts$nums[i] = counts$nums[i] + lengths(DEGs[DEGs$avg_logFC < 0, 1])[[1]]
          counts$nums[i+length(methods)] = counts$nums[i+length(methods)] + lengths(DEGs[DEGs$avg_logFC > 0, 1])[[1]]
        }
      }

      unique_degs[[label]] = unique(all_degs)
      unique_degs$counts[i] = length(unique_degs[[label]])
    }
  })
    
  # adds new variable containing method split
  counts$method_split = factor(str_wrap(unfactor(counts$method),width=12), levels = str_wrap(files$labels,width=12))
  counts$expression = factor(counts$expression, levels = c("Over-Expressed", "Under-Expressed"))

  # plotting summary plot
  p = ggplot(counts, aes(x=method_split, y=nums,fill=expression)) +
    geom_bar(stat="identity", position=position_dodge()) + 
    xlab("Computational Method") + 
    ylab("Number of Differentially Expressed Genes") +
    labs(fill="Differential Expression in\nMethanol-Fixed Samples\nCompared to Fresh Samples")

  # output
  ggsave(output_path[1],p,width=8, height=4)
    
  
  
  # processing for saving as xlsx
  counts$method_split = NULL
  counts["Method"] = counts$method
  
  # adding unique counts to `counts`
  counts["Unique DEGs"] = c(unique_degs$counts, rep(-1, n))
  unique_degs$counts = NULL
  
  # changing format
  counts["Over-Expressed DEGs"] = c(counts$nums[counts$expression == "Over-Expressed"], rep(-1, n))
  counts["Under-Expressed DEGs"] = c(counts$nums[counts$expression == "Under-Expressed"], rep(-1, n))
  
  # removing nums, method and expression columns
  counts["method"] = NULL
  counts["nums"] = NULL
  counts["expression"] = NULL
  counts = counts[1:n, ]
  
  # saving counts to xlsx doc
  wb <- createWorkbook()
  
  s <- createSheet(wb, "Summary")
  addDataFrame(counts, s, row.names = F, col.names = T)
  
  # looping through methods, adding to wb 1 sheet at a time
  labels = gsub(":", "-", names(unique_degs))
  for (i in seq(length(labels))) {
    s <- createSheet(wb, labels[i])
    
    df = data.frame(unique_degs[[i]])
    names(df) = c(paste("UNIQUE.DEGS (", counts[i, 2],")",sep="")) # including number of unique degs in col header
    
    addDataFrame(df, s, row.names = F, col.names = T)
  }
  
  saveWorkbook(wb, output_path[2])
}



################################################################################################
# Concats ARI & NMI into 1 df and saves it
################################################################################################
concat_ari_nmi <- function () {
  # variables
  input_path = config$output_dir
  output_path = paste(config$output_dir, "summary", "ARI_NMI_Summary.xlsx", sep="/")
  methods = config$methods
  n = length(methods)
  
  # reading in files
  files = lapply(methods, function(x) {
    return(read_excel(paste(input_path, x, "clus_contingency_table.xlsx", sep="/"), sheet=1))
  })
  
  stat = c("ARI", "NMI")
  pres = c("", "Fresh", "MeOH")
  df = data.frame("Method" = methods)
  
  # loops through indexes from files
  for (i in seq(2)) {
    for (j in c(2,3)) {
      #  do.call makes a vector of all values from the files list that are on the index [i,j]
      df[paste(stat[i], pres[j], sep=".")] = do.call(c, sapply(files, function(x) return(x[i,j]), USE.NAMES=FALSE))
    }
  }
                                                               
  # saving to xlsx
  wb <- createWorkbook()
  s <- createSheet(wb, "ARI_NMI")
  addDataFrame(df, s, row.names = F, col.names = T)
  addDataFrame(t(data.frame("2" = "* 'No Decontamination' compares cell annotations from the Kidney paper to annotations prior to decontamination",
                            "3" = "* All other methods compare cell annotations prior and post decontamination")),
               s, row.names = F, col.names = F, startRow = (length(df[,1]) + 3))
                                                               
  saveWorkbook(wb, output_path)
                                                               
  return(df)
}


                                                               
################################################################################################
# Concats ARI & NMI into 1 df and saves it
################################################################################################
plot_ari_nmi <- function (df, output_path) {
  # variables
  df = ari_nmi
  output_path = c(paste(config$output_dir, "summary", "ARI_Histogram.png", sep="/"),
                  paste(config$output_dir, "summary", "NMI_Histogram.png", sep="/"))
  n = length(df$Method)

  # 'melting' df so that all values are in 1 col
  new_df = data.frame("method" = rep(df$Method, 2),
                      "stat" = c(rep("ARI", 2*n), rep("NMI", 2*n)),
                      "preservation" = rep(c(rep("Fresh", n), rep("MeOH", n)), 2),
                      "value" = as.numeric(c(df[,2], df[,3], df[,4], df[,5])))
  
  for (i in seq(2)) {
    s = c("ARI", "NMI")[i]
    
    # creating ggplot
    p <- ggplot(new_df[new_df$stat==s,], aes(x=method, y=value, group=preservation)) +
           geom_bar(aes(fill=preservation), stat="identity", position=position_dodge()) + 
           xlab("Computational Method") + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
           scale_y_continuous(name="Value", limits=c(NA, NA))+
           labs(fill="Sample set") + 
           facet_zoom(ylim=c(if (s == "NMI") 0.85 else 0.95,1)) # facet_zoom provides zoomed in section
             
    ggsave(output_path[i], p, height = 6, width = 10)
  }
}


