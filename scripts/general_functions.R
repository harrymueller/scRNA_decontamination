# GENERAL FUNCTIONS
load_libraries <- function () {
  required_libs = c(
    c("Seurat", "readxl", "varhandle", "MASS", "dplyr", "tidyverse"),
    if ("decontaminate" %in% config$process) c("SoupX", "celda") else NULL,
    if ("analyse" %in% config$process) c("ggplot2", "xlsx", "reshape2","plotly","cowplot","patchwork") else NULL
  )
  uninstalled = c()
  
  for (i in required_libs) {
    if (!requireNamespace(i, quietly=T)) {
      # lib not found
      uninstalled = append(uninstalled, i)
    } else {
      # loading lib
      suppressMessages({suppressWarnings({
        library(i, character.only=T)
      })})
    }
  }
  
  if (!is.null(uninstalled)) {
    stop(paste("The following packages are required but not installed:",paste(uninstalled, collapse=", ")))
  } else {
    config = data.frame(quiet=F)
    log_print("All libraries loaded successfully")
  }
}

get_config <- function(args, testing=F) {
  if (!testing)
    config <- config::get(config = (if (length(args) > 0) args[[1]]), file = "config.yml")
  else
    config <- config::get(config="harry")
  
  # TODO check for empty strings in critical places
  
  # splitting "arrays"
  for (i in c("sample_ids", "genes_>=9_cts", "ct_order_dotplots", "pie_plot_cts")) {
    config[[i]] = strsplit(config[[i]], split="|",fixed=T)[[1]]
  }
  return(config)
}

log_print <- function (msg) {
  # TODO output file as well?
  if (!config$quiet)
    print(msg)
}

get_files <- function (config, current_method) {
  files = list(
    "CellRanger" = paste(config$input_dir, "CellRanger", sep="/"),
    "Filtered" = paste(config$input_dir, "Filtered_Feature_Barcode_Matrices", sep="/"),
    "CellBender" = paste(config$input_dir, "CellBender", sep="/"),
    "CellAnnotations" = if (config$recluster == F || (current_method == "none" || current_method == "no_decontamination")) 
                          paste(config$input_dir, config$original_cell_annotations, sep="/") 
                        else paste(config$input_dir, config$new_cell_annotations, sep="/"),
    "GeneSignatures" = paste(config$input_dir, config$gene_signatures, sep="/"),
    "output" = paste(config$output_dir, current_method, sep="/"),
    "OcraRelDir" = paste(config$ocra_dir, current_method, sep="/")
  )
  
  if (substring(current_method,0,5) == "soupx") {
    files$dir = sapply(config$sample_ids, FUN = function(x) {
      paste(files$CellRanger, x, sep="/")
    }, USE.NAMES=F)
    
    files$special = rep(files$GeneSignatures, length(config$sample_ids))
  } else if (substring(current_method,0,7) == "decontx") {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  } else if (current_method == "cellbender") {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$CellBender, "/", x, "/", x, "_filtered.h5", sep="")
    }, USE.NAMES=F)
    
    files$special = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  } else {
    files$dir = sapply(config$sample_ids, function(x) {
      paste(files$CellRanger, x, sep="/")
    }, USE.NAMES=F)
    
    files$special = sapply(config$sample_ids, function(x) {
      paste(files$Filtered,"/", x,".txt", sep="")
    }, USE.NAMES=F)
  }
  
  return(files)
}

get_clusters <- function(path, sample_id, use_new=FALSE) {
  if (!use_new) {
    is_preserved = length(str_split(sample_id, "_",simplify=TRUE))>2
    
    # Annotations for preserved cells are on a separate sheet within the file
    if (!is_preserved) {
      suppressMessages({
        cell_annotations = read_excel(path, sheet=1)
      })
      # library names are different from preserved library names
      cell_annotations = cell_annotations[cell_annotations$Library==sample_id,]
      
    } else {
      suppressMessages({
        cell_annotations = read_excel(path, sheet=2)
      })
      cell_annotations = cell_annotations[cell_annotations$Preservation == "MeOH",c(1:3,5)]
      
      cell_annotations = cell_annotations[cell_annotations$Library==str_split(sample_id, "_",simplify=TRUE)[1],]
    }
    
    names(cell_annotations)[1] = "barcodes"
    names(cell_annotations)[4] = "orig.ident"
    
    cell_annotations$barcodes = sapply(cell_annotations$barcodes, function(x) paste(substring(x, nchar(sample_id)+2), "-1", sep=""))
    names(cell_annotations$orig.ident) = cell_annotations$barcodes
    
    #cell_annotations = subset(cell_annotations, select=-c(1,2,3))
    return(cell_annotations$orig.ident)
  } else {
    annotations = read.table(path, sep="\t",header=F,skip=1,stringsAsFactors=F, as.is=T)
    
    names(annotations) = c("barcode", "celltype")
    
    # creating sample_id variable
    annotations$sample_id = sapply(annotations$barcode, FUN=function(x) {
      return(paste(head(str_split(x,"_",)[[1]],-1),collapse = "_"))
    })
    
    annotations = annotations[which(annotations$sample_id==sample_id),]

    
    # removing sample_id from barcode
    annotations$barcode = sapply(annotations$barcode, FUN=function(x) {
      return(tail(str_split(x, "_")[[1]],1))
    })
    
    ct = annotations$celltype
    names(ct) = annotations$barcode
    return(ct)
  }
}


################################################################################################
### Adds some metadata to the combined seurat object
################################################################################################
adding_metadata <- function(samples.combined) {
  #order_paper <- rev(c("Fib", "MPH", "Podo", "aLOH", "Unknown","CD_PC", "T","NK","DCT","CNT","B","MC","CD_IC","Endo","PT","CD_Trans"))
  Idents(samples.combined)[which(Idents(samples.combined)=="CD_Trans")] = "Unknown"
  order_paper = config$ct_order_dotplots[which(config$ct_order_dotplots %in% levels(samples.combined))]
  
  #samples.combined <- readRDS(paste(output, "samples_integrated_rd.Rda", sep="/"))
  DefaultAssay(samples.combined) <- "RNA"
  
  # adding metadata for `celltype_method`, `celltype`, and changing default ident to `celltype_method`
  samples.combined$celltype_method <- paste(Idents(samples.combined), samples.combined$preservation, sep = "_")
  samples.combined$celltype <- Idents(samples.combined)
  Idents(samples.combined) <- "celltype_method"
  
  samples.combined$celltype.fresh = unlist(lapply(samples.combined$celltype_method, function(x) {
    end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="fresh")
    x <- if (end) str_sub(x,end=-7) else x
  }))
  samples.combined$celltype.meoh = unlist(lapply(samples.combined$celltype_method, function(x) {
    end <- (tail(c(str_split(x, fixed("_"),simplify=T)),n=1)=="MeOH")
    x <- if (end) str_sub(x,end=-6) else x
  }))
  
  return(samples.combined)
}
