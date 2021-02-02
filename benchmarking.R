################################################################################################
# SETUP
################################################################################################
args = commandArgs(trailingOnly=TRUE)

config <- get_config("benchmarking")
config$sample_ids = args[[1]]
current_method = args[[2]]

load_libraries()

# loading functions from separate scripts
source("../improved_scripts/scripts/general_functions.R")
source("../improved_scripts/scripts/decontamination_functions.R")

files=get_files(config, current_method)



################################################################################################
# START
################################################################################################
print(files)
# Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
sample = get_sample(sample_id, files$dir, current_method, files$CellAnnotations, files$special, !(config$recluster == F || (current_method == "none" || current_method == "no_decontamination")))

# ensuring formatting of cell barcodes is the same (across all analyses)
sample$seurat = fix_barcodes(samples[[sample_id]]$seurat)

save_matrices(samples, paste(files$output, "/matrices/",sep=""))

