


### identify transcript origin
### tally transcripts by ct
### express as percentage
### CI?
identify_transcript_origin <- function (samples.combined) {
    # list of genes and origin
    genes = sapply(rownames(samples.combined), function(x) { str_split(x, "-")[[1]][1] })

    # matrix subsets containing only mm || hg genes
    mm_genes = samples.combined@assays$RNA@counts[names(genes)[which(genes == "mm10")],]
    hg_genes = samples.combined@assays$RNA@counts[names(genes)[which(genes == "hg19")],]

    # DF (barcode, ct, endo_counts, non_endo_counts, endo_frac)
    transcripts = data.frame("barcode"=names(Idents(samples.combined)), "celltype"=Idents(samples.combined))
    transcripts$barcode = unfactor(transcripts$barcode)

    # template
    transcripts$native_counts = rep(NaN, length(transcripts$barcode))
    transcripts$non_native_counts = transcripts$native_counts

    # for each subset: for each cell: sum endo transcripts against non-endo transcripts -> calculate fraction of endo
    mm_counts = colSums(as.matrix(mm_genes))
    hg_counts = colSums(as.matrix(hg_genes))

    # set native and non_native counts
    for (i in seq(length(transcripts$barcode))) {
    if (transcripts$celltype[i] == "hg19") {
        transcripts$native_counts[i] = hg_counts[transcripts$barcode[i]]
        transcripts$non_native_counts[i] = mm_counts[transcripts$barcode[i]]
    } else {
        transcripts$native_counts[i] = mm_counts[transcripts$barcode[i]]
        transcripts$non_native_counts[i] = hg_counts[transcripts$barcode[i]]
    }
    }

    transcripts$native_frac = transcripts$native_counts / (transcripts$native_counts + transcripts$non_native_counts)
    transcripts$non_native_frac = 1 - transcripts$native_frac

    return(transcripts)
}

 summarise_counts = function(transcripts) {
    # by ct: average fractions + sd
    # FOCUSES ON NON_NATIVE transcripts ie contamination
    summary_transcripts = data.frame("celltypes" = c("hg19", "mm10", "TOTAL"), "mean"=rep(NaN,3), "std"=rep(NaN,3), "min_ci"=rep(NaN,3),"max_ci"=rep(NaN,3))
    summary_transcripts$celltypes = unfactor(summary_transcripts$celltypes)

    z = qnorm(1-config$alpha/2) # for calculating CI

    # loop through each "ct" (inc. overall)
    for (ct in summary_transcripts$celltypes) {
    # create subsets of transcripts
    if (ct != "TOTAL") {
        transcripts_subset = transcripts[transcripts$celltype == ct,]
        i = if (ct == "hg19") 1 else 2 # index for summary_transcripts
    } else {
        transcripts_subset = transcripts
        i = 3
    }
    
    # mean and standard deviation
    summary_transcripts$mean[i] = mean(transcripts_subset$non_native_frac)
    summary_transcripts$std[i] = sd(transcripts_subset$non_native_frac)
    
    # CI error
    err = z * summary_transcripts$std[i] / sqrt(length(transcripts_subset))
    summary_transcripts$min_ci[i] = summary_transcripts$mean[i] - err
    summary_transcripts$max_ci[i] = summary_transcripts$mean[i] + err
    }

    new = cbind(data.frame("Contaminatio.Fraction.Statistics"=rep("",3)), summary_transcripts)
    names(new) = c("> Contamination Fraction Statistics (as determined by non-endogenous transcripts)", "Cell Types", "Mean", "Std", "Min CI", "Max CI")

    # save to excel: first sheet is summary, then full DF (separated by CT)
    wb <- createWorkbook()
    
    s <- createSheet(wb, "summary")
    addDataFrame(t(new), s, row.names = TRUE, col.names = FALSE)

    for (ct in c("hg19", "mm10")) {
    s = createSheet(wb, ct)
    addDataFrame(transcripts[transcripts$celltype == ct,], s, row.names = FALSE, col.names = TRUE)
    }

    saveWorkbook(wb, "/home/harry/Downloads/test.xlsx")
}

### human gene counts in mouse cells w/ mouse genes x1000 on y && reverse
plot_transcripts <- function(transcripts) {
    # dotplot of transcripts
    p <- ggplot(transcripts, aes(x=native_counts, y=non_native_counts, color=celltype)) + geom_point() +
    theme(text=element_text(size=16, family="TT Times New Roman"))

    ggsave("~/Downloads/counts.png", p, width=10,height=7)
}


