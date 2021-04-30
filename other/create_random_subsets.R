##############################
# USER INPUTS
##############################
N_SUBSETS       = 20
PROP_TO_KEEP    = 0.8

# path to dir containing `genes.tsv` and `barcodes.tsv` 
INPUT_FILE_DIR  = "/data/Perkins/Mouse_Kidney/data/CellRanger/BG3_BG21C/filtered_gene_bc_matrices/mm10"
OUTPUT_FILE_DIR = "/data/Perkins/subset_testing"



##############################
# MAIN SCRIPT
##############################
# repeat for genes and barcodes
for (type in c("genes", "barcodes")) {
    # input
    original = read.table(paste(INPUT_FILE_DIR, "/", type, ".tsv", sep=""), sep="\t")
    original = original[[( if (type == "genes") 2 else 1)]] # only keep relevant column (1 for barcodes && 2 for genes)

    num_to_keep = round(PROP_TO_KEEP * length(original))    # num entries to keep
    output_df = data.frame("V1" = rep("", num_to_keep))     # create template df with template column of correct length

    # loop through N_SUBSETS -> sample -> add to template_df
    for (subset_i in seq(N_SUBSETS)) {
        subset = sample(original, num_to_keep)

        # reorder according to the original order
        subset = original[which(original %in% subset)]  
        output_df[[subset_i]] = subset
    }

    # output to a csv
    write.csv(output_df, paste(OUTPUT_FILE_DIR, "/", type, "_subsets.csv", sep=""),
              row.names = F, col.names = F)
}