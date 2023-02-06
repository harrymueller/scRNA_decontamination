#Libraries
library(Seurat)
library(Matrix)
library(varhandle)
library(stringr)

source("scripts/clustering.R")

# cannabilised from /other/merge_hgmm.R
input_dir = "/data/decont_project/mouse_kidney/input_pt_endo"
#input_dir = "/data/Perkins/Mouse_Kidney/data_endo_pt"
cell_annotations = "/data/decont_project/mouse_kidney/input_pt_endo/new_clus.tsv"
#cell_annotations = "/data/Perkins/Mouse_Kidney/data_endo_pt/new_clus.tsv"
sample_ids = c("BG1_BG20C", "BG3_BG21C", "BG5_BG22C", "BG52_BG20C_MeOH", "BG54_BG21C_MeOH", "BG56_BG22C_MeOH")

cellranger <- function(sample_id, type) {
    # read in data
    dir = paste(input_dir, "CellRanger", sample_id, paste(type, "gene_bc_matrices", sep="_"), "mm10", sep="/")
    sample = Read10X(dir)
    sample = CreateSeuratObject(sample)

    # add clusters
    cell_annotations = get_clusters(cell_annotations, sample_id, FALSE)
    
    print(sample)
    
    # barcodes to remove from the sample
    to_delete = names(cell_annotations[cell_annotations != "PT" & cell_annotations != "Endo"])
    sample = sample[, !colnames(sample) %in% to_delete]
    
    print(sample)
    # output
    writeMM(sample@assays$RNA@counts, paste(dir, "matrix.mtx", sep="/"))

    # barcode file
    write.table(colnames(sample), paste(dir, "barcodes.tsv", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")

    # merge gene
    genes = read.table(paste(dir, "genes.tsv", sep="/"), stringsAsFactor = F)
    genes = genes[which(genes[[2]] %in% rownames(sample)),]
    
    write.table(genes, paste(dir, "genes.tsv", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
}

matrices <- function(sample_id) {
    dir = paste(input_dir, "Filtered_Feature_Barcode_Matrices", paste(sample_id, ".txt", sep=""), sep="/")
    sample = as.matrix(read.csv(dir, header = TRUE, sep = "\t"))

    # add clusters
    cell_annotations = get_clusters(cell_annotations, sample_id, FALSE)

    to_delete = names(cell_annotations[cell_annotations != "PT" & cell_annotations != "Endo"])
    to_delete = sapply(to_delete, function(x) {
        paste(sample_id, strsplit(c(x), "-", fixed=T)[[1]][1], sep = "_")
    })
    names(to_delete) = NULL
    sample = sample[, !colnames(sample) %in% to_delete]

    write.table(as.matrix(sample), file=dir, quote=FALSE, sep="\t")
}

for (s in sample_ids) {
    for (t in c("filtered", "raw")) {
        #cellranger(s, t)
    }
    print(s)
    matrices(s)
}
