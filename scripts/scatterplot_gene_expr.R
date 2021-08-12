gene_expr_scatter_plots = function(prior, post) {
  DefaultAssay(prior) = "RNA"
  
  # for each pres
  for (pres in c("fresh", "MeOH")) {
    # get average expr
    prior_expr = AverageExpression(prior[,prior$preservation == pres], assay = "RNA")$RNA
    post_expr = AverageExpression(post[,post$preservation == pres], assay = "RNA")$RNA
    
    for (ct in names(post_expr)) {
      new = data.frame('X' = log(prior_expr[ct] + 1), 'Y' = log(post_expr[ct] + 1))
      names(new) = c("X", "Y")

      p <- ggplot(new, aes(x=X, y=Y)) +  
        geom_point(alpha = 0.5) + ggtitle(paste("Average gene expression plotted on log(x+1) scales for", pres, "&", ct)) +
        geom_abline(slope=1, intercept=0) + 
        ylab("Post-Decontamination") + xlab("Pre-Decontamination") +
        theme(text=element_text(size=12, family="TT Times New Roman"))
      
      ggsave(paste(files$output, "/plots/gene_expression/", ct, "_", pres, ".png", sep=""), p, width=8, height=8)
    }
  }
}