################################################################################################
# SETUP
################################################################################################
args = commandArgs(trailingOnly=TRUE)

# sample id and current method from arguments
sample_id = args[[1]]
current_method = args[[2]]
i = as.integer(args[[3]])+1

print(paste(rep("#", 25),collapse=""))
print(paste("Starting", current_method, sample_id))

# loading functions from separate scripts
source("scripts/general_functions.R")
source("scripts/decontamination_functions.R")
source("scripts/clustering.R")

# getting config
config <- get_config("benchmarking", "config.yml")
config$sample_ids = sample_id
config$methods=c(current_method)

# fixing config$is_xlsx
config$is_xlsx = sapply(config$methods, function(x) {
  return(config$recluster == F || x == "none" || x == "no_decontamination")
})

# libs
load_libraries()

# files
files=get_files(config, current_method)



################################################################################################
# START
################################################################################################
print("Decontaminating")
# Creates and saves individual R list objects <- previously used `soupx_processing.R` to create the Rda for each sample
sample = get_sample(i, sample_id, current_method)

# ensuring formatting of cell barcodes is the same (across all analyses)
sample$seurat = fix_barcodes(sample$seurat)

save_matrices(list(sample))
print("Completed")
