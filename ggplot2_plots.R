require(ggplot2)

##Code for creating ggplot2 versions of the paper figures.

###############################
#Figure 1 of the supp material#
###############################

##Extract the mutant strains in order
##Change =="m" to "w" for wild type
mut_expr_raw = exprs(yeast.raw[,subset(exp_fac,strain=="m")$data_order])

##Create a ggplot2 friendly data frame
##The data is in the correct order
dd_g = data.frame(m=as.vector(mut_expr_raw), 
                  x = 1:496, 
                  y=rep(1:496, each=496))
dd_g$tps = rep(c(0, 60, 120, 180, 240), each=496^2)
dd_g$replicate = rep(1:3, each=496^2*5)

head(dd_g)
##ggplot2 function
##opts removes axis labels
##scale_ removes tick marks and border
ggplot(dd_g) + geom_raster(aes(x,y,fill=log2(m)), show_guide=FALSE) + 
    facet_grid(replicate ~ tps) + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
         axis.text.y=element_blank(), axis.title.y=element_blank()) +
    scale_x_continuous(expand=c(0, 0), breaks=NULL) + 
    scale_y_continuous(expand=c(0, 0), breaks=NULL)

#####################################
##Figure 2 in supplemental material##
#####################################
#By default, the intensities are only plotted for the pm probes.
Index = unlist(indexProbes(yeast.raw[,1:30], which = "pm"))
intensities = as.vector(intensity(yeast.raw[,1:30])[Index, , drop = FALSE])
no_of_probes = length(intensities)/30

#Construct a data frame for ggplot2
dd_g = data.frame(intensities=log2(intensities), 
        tps = rep(exp_fac$tps, each=no_of_probes),
        replicate = factor(rep(exp_fac$replicate, each=no_of_probes)),
        strain = rep(exp_fac$strain, each=no_of_probes))
dd_g$strain = factor(dd_g$strain, labels=c("Mutant", "Wild type"))

ggplot(dd_g) + 
    geom_density(aes(x=intensities, colour=replicate)) + 
    facet_grid(tps ~ strain) + 
    ylab("Density") + 
    xlab("Log (base 2) intensities") 

#####################################
##Figure 3 in supplemental material##
#####################################
##getMethod("boxplot","AffyBatch")
##getMethod("boxplot","ExpressionSet")
##The boxplots also have range set to 0.
##This corresponds to Inf in ggplot2. 
##Default is 1.5
Index = unlist(indexProbes(yeast.raw, which = "pm"))
intensities = as.vector(intensity(yeast.raw)[Index, , drop = FALSE])
no_of_probes = length(intensities)/30

#Construct a data frame for ggplot2
dd_raw = data.frame(intensities=log2(intensities),
    tps = rep(exp_fac$tps, each=no_of_probes),
    replicate = factor(rep(exp_fac$replicate, each=no_of_probes)),
    strain = rep(exp_fac$strain, each=no_of_probes))
dd_raw$type = "Raw"

no_of_probes = length(as.vector(yeast.matrix))/30
dd_rma = data.frame(intensities=as.vector(yeast.matrix),
                    tps = rep(exp_fac$tps, each=no_of_probes),
                    replicate = factor(rep(exp_fac$replicate, 
                                           each=no_of_probes)),
                    strain = rep(exp_fac$strain, each=no_of_probes))
dd_rma$type = "RMA"


##Boxplots of dd_raw
ggplot(dd_raw) + 
    stat_boxplot(aes(y=intensities, x=replicate), coef=Inf) +
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab("Log (base 2) intensities") + 
    ylim(c(2, 16))

##Boxplots of dd_rma
ggplot(dd_rma) + 
    stat_boxplot(aes(y=intensities, x=replicate), coef=100) + 
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab("Log (base 2) intensities")  + 
    ylim(c(2, 16))

##Boxplots of box raw and rma
ggplot(rbind(dd_raw, dd_rma)) + 
    stat_boxplot(aes(y=intensities, x=replicate,  colour=type), coef=100) + 
    facet_grid(tps ~ strain) +  
    ylab("Density") + 
    xlab("Log (base 2) intensities")


############################################
#Figure 4: timecourse package ##############
############################################
gene_positions = MB.2D$pos.HotellingT2[1:1]
gnames = rownames(yeast.matrix)
gene_probes = gnames[gene_positions]

dd = data.frame(values = as.vector(yeast.matrix[gene_positions,]), 
                strain=exp_fac$strain, 
                tps = exp_fac$tps, 
                replicate = exp_fac$replicate)

ggplot(dd, aes(x=tps, y=values)) + 
    geom_point(aes(colour=strain, shape=factor(replicate)))

plot(yeast.matrix[gene_positions,1:5])






