#Libraries
library(Seurat)
library(Matrix)
source("scripts/hgmm12k/get_data.R")

input_dir = "/media/harry/Data/Perkins/Human_Mouse/data/CellRanger"
output_dir = "/media/harry/Data/Perkins/Human_Mouse/data/CellRanger_merged/raw_gene_bc_matrices"

# INPUT (get data, save as mtx file)
raw = get_raw_hgmm(input_dir, c("hg19", "mm10"))
r = raw@assays$RNA@counts
writeMM(r, paste(output_dir, "matrix.mtx", sep="/"))

# barcode files are identical - only need to copy 1
system2("mv", c(paste(input_dir, "raw_gene_bc_matrices/mm10/barcodes.tsv", sep="/"),
                paste(output_dir, "barcodes.tsv", sep="/")))

# merge genes
system2("cat", c(paste(input_dir, "raw_gene_bc_matrices/mm10/genes.tsv", sep="/"), ">",
                 paste(output_dir, "genes.tsv", sep="/")))
system2("cat", c(paste(input_dir, "raw_gene_bc_matrices/hg19/genes.tsv", sep="/"), ">>", 
                 paste(output_dir, "genes.tsv", sep="/")))

print("Completed.")