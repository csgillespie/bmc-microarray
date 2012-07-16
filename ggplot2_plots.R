library(ggplot2)

##Code for creating ggplot2 versions of the paper figures.

###############################
#Figure 1 of the supp material#
###############################

##Extract the mutant strains in order
##Change =="m" to "w" for wild type
mut_expr_raw = exprs(yeast.raw[,subset(exp_fac,strain=="m")$data_order])

##Create a ggplot2 friendly data frame
##The data is in the correct order
dd_g = data.frame(as.vector(mut_expr_raw), 
                            x = 1:496, 
                            y=rep(1:496, each=496))
dd_g$tps = length(rep(c(0, 60, 120, 180, 240), each=496^2))
dd_g$replicate = length(rep(1:3, each=496^2*5))
                  
##ggplot2 function
##opts removes axis labels
##scale_ removes tick marks and border
ggplot(dd_g) + geom_raster(aes(x,y,fill=log2(m)), show_guide=FALSE) + 
    facet_grid(replicate ~ tps) + 
    opts(axis.text.x=theme_blank(), axis.title.x=theme_blank(), 
         axis.text.y=theme_blank(), axis.title.y=theme_blank()) +
             scale_x_continuous(expand=c(0, 0),breaks=NA) + 
             scale_y_continuous(expand=c(0, 0), breaks=NA)
