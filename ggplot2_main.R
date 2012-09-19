##Figure 1 of main document
##PCA 
dd_pca = exp_fac
dd_pca$tps = factor(dd_pca$tps)
dd_pca$x = yeast.scores[,1]
dd_pca$y = yeast.scores[,2]
ggplot(dd_pca, aes(x, y)) + 
    geom_point(aes(colour=strain, pch=tps))

##Figure 2 of main document
##PCA 
