require(ggplot2)
require(plyr)
require(reshape2)

#############################
##Figure 1 of main document
##PCA 
############################

dd_pca = exp_fac
dd_pca$tps = factor(dd_pca$tps)
dd_pca$x = yeast.scores[,1]
dd_pca$y = yeast.scores[,2]

f1 = ggplot(dd_pca, aes(x, y)) + 
    geom_point(aes(colour=strain, pch=tps)) + 
    xlab("PC 1") + ylab("PC 2") + 
    scale_shape_discrete(name = "Time\n points") +
    scale_colour_discrete(name = "Strain")

#########################################
##Figure 2 of main document
##Time course expression levels for the
##top 9 differently expressed genes
#########################################
top_9_genes = modFordered[1:9]
gene_id = c_df$probe[top_9_genes]
name = c_df$genename[top_9_genes]
exprs_row = yeast.matrix[top_9_genes,]

dd_top_9 = melt(exprs_row) 
dd_top_9$Var2 = as.character(dd_top_9$Var2)
##Calculate the means values over replicates
##Also calculate ymin and ymax in case we want error bars
dd_top_9 = ddply(dd_top_9, .(Var1, Var2), summarise, 
        m = mean(value), ymax = max(value), ymin=min(value),
        name = name[gene_id==Var1])

##Take the first character for strain
dd_top_9$Strain = substr(dd_top_9$Var2, 1, 1)
dd_top_9$tps = substr(dd_top_9$Var2, 2, nchar(dd_top_9$Var2))
dd_top_9$tps = as.numeric(dd_top_9$tps)
dd_top_9$name = factor(dd_top_9$name, levels=name)

f2 = ggplot(dd_top_9, aes(tps)) + 
    geom_line(aes(y=m, colour=Strain)) +
    geom_point(aes(y=m, colour=Strain)) +
    facet_wrap(~name, ncol=3) + 
    xlab("Time") + 
    ylab("Expression level")

######################################################
##Figure 3 of volcano plot
##Highlight genes that have an absolute 
##fold change > 2 and a p-value < Bonferroni cut-off
######################################################

dd_vol = topTable(eb, coef=1, number=1000000, sort.by="logFC")
no_of_genes = nrow(dd_vol)
dd_vol$threshold = as.factor(abs(gene_list$logFC) > log2(2) 
                                & gene_list$P.Value < 0.05/no_of_genes)

## Construct the plot object
f3 = ggplot(dd_vol, aes(logFC, -log10(P.Value), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) +
    theme(legend.position = "none") +
    geom_hline(yintercept=-log10(0.05/5900), alpha=0.3) + 
    geom_vline(xintercept=log2(2), alpha=0.3) + 
    geom_vline(xintercept=-log2(2), alpha=0.3) + 
    xlab(expression(paste(log[2], " fold change"))) + 
    ylab(expression(paste(-log[10], " p-value"))) +
    ylim(c(0, 30))



