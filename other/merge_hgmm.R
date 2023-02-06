#Libraries
library(Seurat)
library(Matrix)
library(varhandle)

source("/data/Perkins/develop/scripts/hgmm12k/get_data.R")

input_dir = "/data/Perkins/hgmm12k/data/CellRanger"

for (type in c("raw", "filtered")) {
    output_dir = paste("/data/Perkins/hgmm12k/data/CellRanger_merged2/", type, "_gene_bc_matrices", sep="")

    # INPUT (get data, save as mtx file)
    if (type == "filtered")
        raw = get_filtered_hgmm(input_dir, "/data/Perkins/hgmm12k/data/gem_classification.csv", c("hg19", "mm10"))
    else
        raw = get_raw_hgmm(input_dir, c("hg19", "mm10"))

    writeMM(raw@assays$RNA@counts, paste(output_dir, "matrix.mtx", sep="/"))

    # barcode files are identical - only need to copy 1
    write.table(colnames(raw), paste(output_dir, "barcodes.tsv", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")

    # merge genes
    genes = rbind(read.table(paste(input_dir, "/", type, "_gene_bc_matrices/mm10/genes.tsv", sep=""), stringsAsFactor = F),
                  read.table(paste(input_dir, "/", type, "_gene_bc_matrices/hg19/genes.tsv", sep=""), stringsAsFactor = F))

    if (type == "filtered")
        genes = genes[which(genes[[2]] %in% rownames(raw)),]
    write.table(genes, paste(output_dir, "genes.tsv", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")

    #system2("cat", c(paste(input_dir, "/", type, "_gene_bc_matrices/hg19/genes.tsv", sep=""), ">", 
    #                paste(output_dir, "genes.tsv", sep="/")))
    #system2("cat", c(paste(input_dir, "/", type, "_gene_bc_matrices/mm10/genes.tsv", sep=""), ">>",
    #               paste(output_dir, "genes.tsv", sep="/")))

}
print("Completed.")
