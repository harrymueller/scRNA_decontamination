################################################################################################
# SETUP
################################################################################################
args = commandArgs(trailingOnly=TRUE)

print(paste("Starting",args[[2]],args[[1]]))


# loading functions from separate scripts
source("scripts/general_functions.R")
source("scripts/decontamination_functions.R")

config <- get_config("benchmarking")
config$sample_ids = args[[1]]
sample_id = args[[1]]
current_method = args[[2]]

load_libraries()


files=get_files(config, current_method)



################################################################################################
# START
################################################################################################
print("Decontaminating")
# Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
sample = get_sample(sample_id, files$dir, current_method, files$CellAnnotations, files$special, !(config$recluster == F || (current_method == "none" || current_method == "no_decontamination")))

# ensuring formatting of cell barcodes is the same (across all analyses)
sample$seurat = fix_barcodes(sample$seurat)

save_matrices(list(sample), paste(files$output, "/matrices/",sep=""))
print("Completed")
