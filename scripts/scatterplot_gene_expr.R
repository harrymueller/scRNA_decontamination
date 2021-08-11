scatterplot_gene_expr_between_samples <- function (obj1, obj2, labs) {
    # by preservation
    for (pres in c("fresh", "MeOH")) {
        # average gene expr
        avg1 = average_expr(obj1, pres)
        avg2 = average_expr(obj2, pres)
        
        avg_expr = data.frame("x" = avg1, "y" = avg2[names(avg1)])
        # scatterplot
        p <- ggplot(avg_expr, aes(x=log(x+1), y=log(y+1))) +  
                  geom_point(alpha = 0.5) + ggtitle(paste("Gene Expression Comparison; log(x+1) scales", pres)) +
                  geom_abline(slope=1, intercept=0) + 
                  ylab(labs[2]) + xlab(labs[1]) +
                  scale_color_manual(values=c("#27AE60", "#8E44AD", "#E67E22","#3498DB"), name="Species and Prior- or\nPost-Decontamination") +
                  theme(text=element_text(size=16, family="TT Times New Roman"))
    
        ggsave(paste(files$output, "/plots/gene_expression_", pres, ".png", sep=""), p, width=8, height=8)
    }
        
}


average_expr <- function (seurat, pres) {
    # returns a named vector of the gene expression averages for the preservation type given
    return(rowMeans(as.matrix(seurat[, seurat$preservation == pres]@assays$RNA@counts)))
}