# Benchmarking algorithms for removing ambient RNA signals from single-cell RNA-seq data

### Authors
Harrison A Mueller<sup>1</sup>, Alistair R R Forrest<sup>1</sup>, Elena Denisenko<sup>1</sup>

*1. Harry Perkins Institute of Medical Research, QEII Medical Centre and Centre for Medical Research, The University of Western Australia, Nedlands, Perth, WA 6009, Australia*

##
This repository stores the scripts used throughout this project.

## Abstract
**Background:** Single-cell RNA sequencing technologies have advanced significantly over the last few years, and are becoming increasingly important in medical research. It has been shown that contamination is present in the form of ambient RNA in the cell buffer. Several methods have been developed to remove this contamination. However, there has been no comparison of the efficacy of these tools.

**Results:** In this paper, four decontamination algorithms (SoupX, DecontX, FastCAR, and CellBender) were applied to two datasets, some with varying parameters. The first dataset was a human-mouse mixed dataset, in which it was trivial to identify removal of inter-species contamination. CellBender was able to remove most of the contamination from the majority of the cells, and both SoupX and DecontX performed well on cells with a large starting contamination. FastCAR did not perform well. The second dataset was a mouse kidney dataset, which was a more accurate representation of real-world datasets. Several analyses were performed, including investigating the changes to cell type annotations, and differential gene expression analyses. Few annotations changed following re-annotation after decontamination, with the exception of FastCAR. SoupX and CellBender were able to remove a large portion of the genes identified as significantly contributing to the contamination. DecontX was able to reduce some of the contamination, and FastCAR was able to reduce some contaminating transcripts heavily, but did not alter the expression of others.

**Conclusions:** CellBender was the most effective algorithm on data that is not heavily contaminated, and SoupX performs well on samples with a large amount of contamination. 

**Keywords:** Single-cell transcriptomics, scRNA-seq, ambient RNA contamination, SoupX, DecontX, FastCAR, CellBender
