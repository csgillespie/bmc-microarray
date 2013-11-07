require(ggplot2)

##Code for creating ggplot2 versions of
##the figures in the supplemental material

#########################################
## Figure 1 of the supp material        #
## Image plots of the probe intensities #
#########################################
##Extract the mutant strains in order
##Change == "m" to "w" for wild type
mut_expr_raw = exprs(yeast.raw[,subset(exp_fac,strain=="m")$data_order])

##Create a ggplot2 friendly data frame
##The data is in the correct order
dd_int = data.frame(m=as.vector(mut_expr_raw), 
                  x = 1:496, 
                  y=rep(1:496, each=496))
dd_int$tps = rep(c(0, 60, 120, 180, 240), each=496^2)
dd_int$replicate = rep(1:3, each=496^2*5)

##theme removes axis labels
##scale_ removes tick marks and border
s1 = ggplot(dd_int) + 
    geom_raster(aes(x,y,fill=log2(m)), show_guide=FALSE) + 
    facet_grid(replicate ~ tps) + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
         axis.text.y=element_blank(), axis.title.y=element_blank()) +
    scale_x_continuous(expand=c(0, 0), breaks=NULL) + 
    scale_y_continuous(expand=c(0, 0), breaks=NULL)

##########################################
##Figure 2 in supplemental material     ##
##Density plots of the log2 intensities ##
##########################################
#By default, the intensities are only plotted for the pm probes.
index = unlist(indexProbes(yeast.raw[,1:30], which = "pm"))
intensities = as.vector(intensity(yeast.raw[,1:30])[index, , drop = FALSE])
no_of_probes = length(intensities)/30

#Construct a data frame for ggplot2
dd_den = data.frame(intensities=log2(intensities), 
        tps = rep(exp_fac$tps, each=no_of_probes),
        replicate = factor(rep(exp_fac$replicate, each=no_of_probes)),
        strain = rep(exp_fac$strain, each=no_of_probes))
dd_den$strain = factor(dd_den$strain, labels=c("Mutant", "Wild type"))

s2 = ggplot(dd_den) + 
    geom_density(aes(x=intensities, colour=replicate)) + 
    facet_grid(tps ~ strain) + 
    ylab("Density") + 
    xlab(expression(paste(log[2], " intensities"))) + 
    theme_bw()
s2

#################################################
##Figure 3 in supplemental material             #
##Boxplot of the raw and normalised intensities #
#################################################
##getMethod("boxplot","AffyBatch")
##getMethod("boxplot","ExpressionSet")

##The boxplots also have range set to 0.
##This corresponds to Inf in ggplot2. 
##Default is 1.5
index = unlist(indexProbes(yeast.raw, which = "pm"))
intensities = as.vector(intensity(yeast.raw)[index, , drop = FALSE])
no_of_probes = length(intensities)/30

#Construct a data frame for ggplot2
dd_box_raw = data.frame(intensities=log2(intensities),
    tps = rep(exp_fac$tps, each=no_of_probes),
    replicate = factor(rep(exp_fac$replicate, each=no_of_probes)),
    strain = rep(exp_fac$strain, each=no_of_probes))
dd_box_raw$type = "Raw"

no_of_probes = length(as.vector(yeast.matrix))/30
dd_box_rma = data.frame(intensities=as.vector(yeast.matrix),
                    tps = rep(exp_fac$tps, each=no_of_probes),
                    replicate = factor(rep(exp_fac$replicate, 
                                           each=no_of_probes)),
                    strain = rep(exp_fac$strain, each=no_of_probes))
dd_box_rma$type = "RMA"

##Boxplots of Raw intensities
s3a = ggplot(dd_box_raw) + 
    stat_boxplot(aes(replicate, intensities), coef=Inf) +
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab(expression(paste(log[2], " intensities"))) + 
    coord_cartesian(ylim=c(2, 16))

##Boxplots of transformed intensities
s3b = ggplot(dd_box_rma) + 
    stat_boxplot(aes(replicate, intensities), coef=Inf) + 
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab(expression(paste(log[2], " intensities"))) + 
    coord_cartesian(ylim=c(2, 16))

##Boxplots of box raw and rma
s3_joint = ggplot(rbind(dd_box_raw, dd_box_rma)) + 
    stat_boxplot(aes(replicate, intensities, colour=type), coef=100) + 
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab(expression(paste(log[2], " intensities"))) + 
    scale_colour_discrete(name = "Type")
    

############################################
##Figure 4: timecourse package            ##
##Top differentially expressed gene       ##
############################################
n = 9
gene_positions = MB.2D$pos.HotellingT2[1:n]
gnames = rownames(yeast.matrix)
gene_probes = gnames[gene_positions]
cn = colnames(yeast.matrix)

head(yeast.matrix)
dd_tc = data.frame(values = as.vector(yeast.matrix[gene_positions,]))
dd_tc$Strain = rep(substr(cn, 1, 1), each=n)
dd_tc$tps = rep(substr(cn, 2, nchar(cn)), each=n)
dd_tc$tps = as.numeric(dd_tc$tps)
dd_tc$Replicate = factor(rep(c(1:3, 1:3), each=5*n))
dd_tc$names = gnames[gene_positions]
dd_tc$names = factor(dd_tc$names, labels = gnames[gene_positions])

ggplot(dd_tc, aes(tps, values)) + 
    geom_point(aes(colour=Strain, shape=Replicate)) + 
    facet_wrap(~names, ncol=3) + 
    xlab("Time") + 
    ylab("Expression level")


