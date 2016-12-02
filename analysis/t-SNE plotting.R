library(Rtsne)
library(ggplot2)

set.seed(47)

tsne_out = Rtsne(as.matrix(expr.retina.var.stab), perplexity = 5, check_duplicates = F, theta = 0.0)
tsne_plot = data.frame(tsne_out$Y)

tsne_plot %>%
  ggplot(aes(x=X1,y=X2,colour=factor(meta.retina$study_accession))) + 
  geom_point(size=4) + ggtitle(paste0("t-sne. Perplexity = ", 5))