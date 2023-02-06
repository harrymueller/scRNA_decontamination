##############################
# Creates random subsets of a given dataset
##############################



##############################
# USER INPUTS
##############################
N_SUBSETS       = 2
PROP_TO_KEEP    = 0.8

# path to dir containing `genes.tsv` and `barcodes.tsv` 
INPUT_FILE_DIR  = c("/data/Perkins/hgmm12k/data/CellRanger_merged/raw_gene_bc_matrices",
                    "/data/Perkins/hgmm12k/data/CellRanger_merged/filtered_gene_bc_matrices")
OUTPUT_FILE_DIR = "/data/Perkins/stability_testing2"

GENERATE_NEW_SUBSETS = T



##############################
# LIBRARIES
##############################
library(varhandle)
library(Seurat)
library(Matrix)



##############################
# Produces random subsets of the input data
    # returns a list containing two DF for gene and barcode subsets
    # also saves the gene and barcode subsets to a CSV
##############################
produce_random_subsets <- function (type) {
    output_df = list()

    # input
    if (type == "genes")
        original = unique(GENES_FULL[[2]])
    else
        original = unique(BARCODES_FULL[[1]])

    num_to_keep = round(PROP_TO_KEEP * length(original))    # num entries to keep
    output_df[[type]] = data.frame("V1" = rep("", num_to_keep))     # create template df with template column of correct length

    # loop through N_SUBSETS -> sample -> add to template_df
    for (subset_i in seq(N_SUBSETS)) {
        subset = sample(original, num_to_keep)

        # reorder according to the original order
        subset = original[which(original %in% subset)]  
        output_df[[type]][[subset_i]] = subset

        if (type == "genes")
            output_df[["barcodes"]][[subset_i]] = BARCODES_FULL$V1
        else
            output_df[["genes"]][[subset_i]] = GENES_FULL$V2
    }

    # output to a csv
    write.csv(output_df[[type]], paste(OUTPUT_FILE_DIR, "/", type, "_subsets.csv", sep=""),
            row.names = F)
    
    

    return(output_df)
}


##############################
# Reads in the subsets saved as CSV
##############################
read_random_subsets <- function (type) {
    subset      = read.csv(paste(OUTPUT_FILE_DIR, "/", type, "_subsets.csv", sep=""),
                           stringsAsFactors = FALSE)

    if (type == "genes")
        return(list("genes" = subset, "barcodes" = BARCODES_FULL))
    else
        return(list("genes" = GENES_FULL, "barcodes" = subset))
}



##############################
# Given an index -> subsets all the files and saves the files
##############################
save_subset <- function(index, type, data_index) {
    # make folder
    output_dir = paste(OUTPUT_FILE_DIR, "/data_", type, "/", index, sep="")
    if (!dir.exists(output_dir))
        dir.create(output_dir)

    if (data_index == 1)
        output_dir = paste(output_dir, "raw_gene_bc_matrices", sep="/")
    else
        output_dir = paste(output_dir, "filtered_gene_bc_matrices", sep="/")

    if (!dir.exists(output_dir))
            dir.create(output_dir)

    # subset and save genes
    genes_subset = GENES_FULL[GENES_FULL$V2 %in% subsets$genes[[index]],]
    write.table(genes_subset, paste(output_dir, "genes.tsv", sep="/"), 
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t") 
    
    # subset and save barcodes
    barcodes_subsets = BARCODES_FULL[BARCODES_FULL$V1 %in% subsets$barcodes[[index]],]
    write.table(barcodes_subsets, paste(output_dir, "barcodes.tsv", sep="/"))    

    matrix_subset = MATRIX_FULL[rownames(MATRIX_FULL) %in% subsets$genes[[index]],
                                colnames(MATRIX_FULL) %in% subsets$barcodes[[index]]]
    writeMM(matrix_subset, paste(output_dir, "matrix.mtx", sep="/"))
}



##############################
# MAIN EXECUTION
##############################
for (data_index in c(1,2)) {
    if (data_index == 1)
        print("Raw")
    else
        print("Filtered")
        
    # load in full data
    BARCODES_FULL = read.table(paste(INPUT_FILE_DIR[[data_index]], "/", "barcodes", ".tsv", sep=""), sep="\t")
    GENES_FULL    = read.table(paste(INPUT_FILE_DIR[[data_index]], "/", "genes", ".tsv", sep=""), sep="\t")
    MATRIX_FULL   = Read10X(INPUT_FILE_DIR[[data_index]])

    # unfactor
    BARCODES_FULL[[1]] = unfactor(BARCODES_FULL[[1]])
    GENES_FULL[[1]] = unfactor(GENES_FULL[[1]])
    GENES_FULL[[2]] = unfactor(GENES_FULL[[2]])

    # generate subsets
    for (type in c("genes", "barcodes")) {
        # ensure data folder is present
        if (!dir.exists(paste(OUTPUT_FILE_DIR, "/data_", type, sep="")))
            dir.create(paste(OUTPUT_FILE_DIR, "/data_", type, sep=""))

        if (GENERATE_NEW_SUBSETS & data_index == 1) {
            print("Generating subsets")
            subsets = produce_random_subsets(type)
        } else {
            print("Loading subsets")
            subsets = read_random_subsets(type)
        }

        # save subsets
        print("Saving subsets")
        for (i in seq(N_SUBSETS))
            save_subset(i, type, data_index)

        print("All subsets saved")
    }
}
