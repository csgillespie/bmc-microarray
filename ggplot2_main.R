library(ggplot2)
library(plyr)
library(reshape2)
##Figure 1 of main document
##PCA 
dd_pca = exp_fac
dd_pca$tps = factor(dd_pca$tps)
dd_pca$x = yeast.scores[,1]
dd_pca$y = yeast.scores[,2]
ggplot(dd_pca, aes(x, y)) + 
    geom_point(aes(colour=strain, pch=tps)) + 
    xlab("PC 1") + ylab("PC 2")

##Figure 2 of main document
##PCA 
top_9_genes = modFordered[1:9]
gene_id = c_df$probe[top_9_genes]
name = c_df$genename[top_9_genes]
exprs.row = yeast.matrix[top_9_genes,]

dd = melt(exprs.row) 
dd$Var2 = as.character(dd$Var2)
dd = ddply(dd, .(Var1, Var2), summarise, 
        m = mean(value), ymax = max(value), ymin=min(value),
        name = name[gene_id==Var1])
dd$Strain = substr(dd$Var2, 1, 1)
dd$tps = substr(dd$Var2, 2, nchar(dd$Var2))
dd$tps = as.numeric(dd$tps)
dd$name =factor(dd$name, levels=name)

ggplot(dd, aes(tps)) + 
    geom_line(aes(y=m, colour=Strain)) +
    geom_point(aes(y=m, colour=Strain)) +
    facet_wrap(~name, ncol=3) + 
    xlab("Time") + ylab("Expression level")

##Figure 3 of volcano plot

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list = topTable(eb, coef=1, number=1000000, sort.by="logFC")
no_of_genes = nrow(gene_list)
gene_list$threshold = as.factor(abs(gene_list$logFC) > log2(2) 
                                & gene_list$P.Value < 0.05/no_of_genes)

##Construct the plot object
ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) +
    opts(legend.position = "none") +
    geom_hline(yintercept=-log10(0.05/5900), alpha=0.3) + 
    geom_vline(xintercept=log2(2), alpha=0.3) + 
    geom_vline(xintercept=-log2(2), alpha=0.3) + 
    xlab(expression(paste(log[2], " fold change"))) + 
    ylab(expression(paste(-log[10], " p-value"))) +
    ylim(c(0, 30))
    
#Figure 4
#Heat map

