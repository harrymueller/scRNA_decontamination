##############################
# USER INPUTS
##############################
N_SUBSETS       = 20
PROP_TO_KEEP    = 0.8

# path to dir containing `genes.tsv` and `barcodes.tsv` 
INPUT_FILE_DIR  = "/data/Perkins/Mouse_Kidney/data/CellRanger/BG3_BG21C/filtered_gene_bc_matrices/mm10"
OUTPUT_FILE_DIR = "/data/Perkins/subset_testing"

GENERATE_NEW_SUBSETS = T



##############################
# LIBRARIES
##############################
library(varhandle)
library(Seurat)
library(Matrix)



##############################
# FUNCTIONS
##############################
produce_random_subsets <- function () {
    output_df = list()

    # repeat for genes and barcodes
    for (type in c("genes", "barcodes")) {
        # input
        if (type == "genes")
            original = GENES_FULL[[2]]
        else
            original = BARCODES_FULL[[1]]

        num_to_keep = round(PROP_TO_KEEP * length(original))    # num entries to keep
        output_df[[type]] = data.frame("V1" = rep("", num_to_keep))     # create template df with template column of correct length

        # loop through N_SUBSETS -> sample -> add to template_df
        for (subset_i in seq(N_SUBSETS)) {
            subset = sample(original, num_to_keep)

            # reorder according to the original order
            subset = original[which(original %in% subset)]  
            output_df[[type]][[subset_i]] = subset
        }

        # output to a csv
        write.csv(output_df[[type]], paste(OUTPUT_FILE_DIR, "/", type, "_subsets.csv", sep=""),
                row.names = F)
    }

    return(output_df)
}



# using another function to read in subsets so you can maintain the same random subset
read_random_subsets <- function () {
    genes       = read.csv(paste(OUTPUT_FILE_DIR, "/", "genes_subsets.csv", sep=""))
    barcodes    = read.csv(paste(OUTPUT_FILE_DIR, "/", "barcodes_subsets.csv", sep=""))

    for (col in colnames(genes)) {
        genes[col] = unfactor(genes[col])
        barcodes[col] = unfactor(barcodes[col])
    }

    return(list("genes" = genes, "barcodes" = barcodes))
}



# given an index - subsets genes.tsv, barcodes.tsv, matrix.mtx and saves
save_subset <- function(index) {
    output_dir = paste(OUTPUT_FILE_DIR, "data", index, sep="/")
    # make folder
    if (!dir.exists(output_dir))
        dir.create(output_dir)

    # subset and save genes
    genes_subset = GENES_FULL[GENES_FULL$V2 %in% subsets$genes[[index]],]
    write.table(genes_subset, paste(output_dir, "genes.tsv", sep="/"), 
                col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t") 

    # subset and save barcodes
    barcodes_subsets = BARCODES_FULL[BARCODES_FULL$V1 %in% subsets$barcodes[[index]],]
    write.table(barcodes_subsets, paste(output_input
    matrix_subset = MATRIX_FULL[rownames(MATRIX_FULL) %in% subsets$genes[[index]],
                                colnames(MATRIX_FULL) %in% subsets$barcodes[[index]]]
    writeMM(matrix_subset, paste(output_dir, "matrix.mtx", sep="/"))
}



##############################
# MAIN SCRIPT
##############################
# load in full data
BARCODES_FULL = read.table(paste(INPUT_FILE_DIR, "/", "barcodes", ".tsv", sep=""), sep="\t")
GENES_FULL    = read.table(paste(INPUT_FILE_DIR, "/", "genes", ".tsv", sep=""), sep="\t")
MATRIX_FULL   = Read10X(INPUT_FILE_DIR)

# unfactor
BARCODES_FULL[[1]] = unfactor(BARCODES_FULL[[1]])
GENES_FULL[[1]] = unfactor(GENES_FULL[[1]])
GENES_FULL[[2]] = unfactor(GENES_FULL[[2]])

# generate subsets
if (GENERATE_NEW_SUBSETS) {
    print("Generating subsets")
    subsets = produce_random_subsets()
} else {
    print("Loading subsets")
    subsets = read_random_subsets()
}

# save subsets
print("Saving subsets")
for (i in seq(N_SUBSETS))
    save_subset(i)

print("All subsets saved")