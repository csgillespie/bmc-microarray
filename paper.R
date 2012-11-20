
## @knitr unnamed-chunk-1
opts_chunk$set(prompt = FALSE,fig.path='graphics/', fig.align='center',fig.lp="")
#knit_hooks$set()
options(width=75)
knit_hooks$set(
    par.nice = function(before, options, envir) {
        if (before && !is.null(options$par.mfrow)){
            par(mar=c(3,3,2,1),  mgp=c(2,0.4,0), tck=-.01,
                cex.lab=.95,cex.axis=0.9, las=1, mfrow=options$par.mfrow)
            
        } else if (before) {
            par(mar=c(3,3,2,1),  mgp=c(2,0.4,0), tck=-.01,
                cex.lab=.95,cex.axis=0.9, las=1)
        }
    }, 
    crop=hook_pdfcrop)
read_chunk('ExtractIDs.R')
read_chunk('RemoveProbes.R')


## @knitr unnamed-chunk-2
url = 'http://bioconductor.org/biocLite.R'
source(url)
biocLite()


## @knitr unnamed-chunk-3
##From Bioconductor
biocLite(c('ArrayExpress', 'Mfuzz', 'timecourse', 'yeast2.db', 
    'yeast2probe', 'yeast2cdf'))
##From cran
install.packages(c('GeneNet', 'gplots'))


## @knitr unnamed-chunk-4
update.packages(repos=biocinstallRepos())


## @knitr unnamed-chunk-5
library(ArrayExpress)
yeast.raw = ArrayExpress('E-MEXP-1551')


## @knitr unnamed-chunk-6
ph = yeast.raw@phenoData
exp_fac = data.frame(data_order = 1:30, 
    strain = ph@data$Factor.Value..GENOTYPE.,  
    replicates = ph@data$Factor.Value..INDIVIDUAL.,
    tps = ph@data$Factor.Value..TIME.)
levels(exp_fac$strain) = c('m', 'w')
exp_fac = with(exp_fac, exp_fac[order(strain, replicates, tps), ])
exp_fac$replicate = rep(1:3, each=5, 2)


## @knitr unnamed-chunk-7
#Read in the mask file
s_cer = read.table('s_cerevisiae.msk', skip=2, 
stringsAsFactors=FALSE)
source('ExtractIDs.R')
c_df = ExtractIDs(s_cer[ ,1])


## @knitr unnamed-chunk-8
#Get the raw dataset for S. cerevisiae only
library(affy)
library(yeast2probe)
source('RemoveProbes.R')
cleancdf = cleancdfname(yeast.raw@cdfName)
RemoveProbes(s_cer[, 1], cleancdf, 'yeast2probe') 


## @knitr s_fig1
for(i in 1:5) {
    plot_title = paste('Strain: ', exp_fac$strain[i], 
        'Time: ', exp_fac$tps[i])
        d = exp_fac$data_order[i]
    image(yeast.raw[,d], main=plot_title) 
}


## @knitr s_fig2
d = exp_fac$data_order[1:5]
hist(yeast.raw[,d], lwd=2, ylab="Density", xlab="Log (base 2) intensities") 


## @knitr unnamed-chunk-9
yeast.rma = rma(yeast.raw)
yeast.matrix = exprs(yeast.rma)[, exp_fac$data_order]
colnames(yeast.matrix) = paste0(exp_fac$strain, exp_fac$tps)
exp_fac$data_order = 1:30


## @knitr s_fig3
library(affyPLM)
#Raw data intensities
boxplot(yeast.raw, col='red', main='', ylim=c(2,16))
#Normalised intensities
boxplot(yeast.rma, col='blue', ylim=c(2,16))


## @knitr unnamed-chunk-10
yeast.PC = prcomp(t(yeast.matrix))
yeast.scores = predict(yeast.PC)


## @knitr F1
#Plot of the first two principal components
plot(yeast.scores[,1], yeast.scores[,2], xlab='PC 1', ylab='PC 2', 
    pch=rep(1:5, 6), col=as.numeric(exp_fac$strain))
legend("bottomleft", pch=1:5, cex=0.6, 
    c('t 0', 't 60', 't 120', 't 180', 't 240'))


## @knitr unnamed-chunk-11
library(timecourse)
size = matrix(3, nrow = 5900, ncol = 2)


## @knitr unnamed-chunk-12
c.grp = as.character(exp_fac$strain)
t.grp = as.numeric(exp_fac$tps)
r.grp = as.character(exp_fac$replicate)
MB.2D = mb.long(yeast.matrix, times = 5, method = '2', 
    reps = size, condition.grp = c.grp, time.grp = t.grp, 
    rep.grp = r.grp)


## @knitr unnamed-chunk-13
gene_positions = MB.2D$pos.HotellingT2[1:100]
gnames = rownames(yeast.matrix)
gene_probes = gnames[gene_positions]


## @knitr s_fig4
plotProfile(MB.2D, ranking=1, gnames=rownames(yeast.matrix))


## @knitr unnamed-chunk-14
library(limma)
expt_structure = factor(colnames(yeast.matrix))

#Construct the design matrix
X = model.matrix(~0 + expt_structure)
colnames(X) =  c('m0', 'm60', 'm120', 'm180', 'm240', 
    'w0', 'w60', 'w120', 'w180', 'w240')


## @knitr unnamed-chunk-15
lm.fit = lmFit(yeast.matrix, X) 


## @knitr unnamed-chunk-16
mc = makeContrasts('m60-w60', 'm120-w120', 'm180-w180', 'm240-w240', levels=X)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)


## @knitr unnamed-chunk-17
#see help(toptable) for more options
toptable(eb, sort.by='logFC') 


## @knitr unnamed-chunk-18
topTableF(eb)


## @knitr unnamed-chunk-19
modFpvalue = eb$F.p.value
##Change 'bonferroni' to 'fdr' to use the 
##false discovery rate as a cut-off
indx = p.adjust(modFpvalue, method='bonferroni') < 0.05
sig = modFpvalue[indx]

#No. of sig. differential expressed genes
nsiggenes = length(sig) 
results = decideTests(eb, method='nestedF')

modF = eb$F
modFordered = order(modF, decreasing = TRUE)

#Retrieve the significant probes and genes
c_rank_probe = c_df$probe[modFordered[1:nsiggenes]]
c_rank_genename = c_df$genename[modFordered[1:nsiggenes]]

#Create a list and write to a file
updown = results[modFordered[1:nsiggenes],]
write.table(cbind(c_rank_probe, c_rank_genename, updown), 
    file='updown.csv', sep=',', 
    row.names=FALSE, col.names=FALSE)


## @knitr F2

for(i in 0:8){
    indx = rank(modF) == nrow(yeast.matrix)-i
    
    id = c_df$probe[indx]
    name = c_df$genename[indx]
    exprs.row = yeast.matrix[indx,]
    genetitle = paste(sprintf('%.30s', id), sprintf('%.30s', name), 
                      'Rank =', i+1) 

    plot(0, pch=NA, xlim=range(0, 240), ylim=range(exprs.row), 
         ylab='Expression', xlab='Time', main=genetitle,
         cex.main=0.7)
    
    for(j in 1:6){
        pch_value = as.character(exp_fac$strain[5*j])
        points(c(0, 60, 120, 180, 240), exprs.row[(5*j-4):(5*j)], 
               type='b', pch=pch_value)
    }  
} 


## @knitr unnamed-chunk-20
#Rank of Probesets, also output gene names
par(mfrow=c(3, 3), ask=TRUE)
for(i in 0:99){
    indx = rank(modF) == nrow(yeast.matrix) - i
    id = c_df$probe[indx]
    name = c_df$genename[indx]
    exprs.row = yeast.matrix[indx,]
    genetitle = paste(sprintf('%.30s', id), sprintf('%.30s', name), 
                      'Rank =', i+1) 
    
    plot(0, pch=NA, xlim=range(0, 240), ylim=range(exprs.row), 
         ylab='Expression', xlab='Time', main=genetitle)
    
    for(j in 1:6){
        pch_value = as.character(exp_fac$strain[5*j])
        points(c(0, 60, 120, 180, 240), exprs.row[(5*j-4):(5*j)], 
               type='b', pch=pch_value)
    }  
} 


## @knitr unnamed-chunk-21
N = 100
gene_positions = MB.2D$pos.HotellingT2[1:N]
tc_top_probes = gnames[gene_positions]
lm_top_probes = c_df$probe[modFordered[1:N]]
length(intersect(tc_top_probes, lm_top_probes))


## @knitr unnamed-chunk-22
#Obtain the maximum fold change but keep the sign 
maxfoldchange = function(foldchange){
    foldchange[which.max(abs(foldchange))]
}
    
difference = apply(eb$coeff, 1, maxfoldchange)
pvalue = eb$F.p.value
lodd = -log10(pvalue)

#hfc: high fold-change
nd = (abs(difference)>log(2, 2))
ordered_hfc = order(abs(difference), decreasing=TRUE)
hfc = ordered_hfc[1:length(difference[nd])]

np = p.adjust(pvalue, method='bonferroni') < 0.05

#lpv: low p-value (large F-value)
ordered_lpv = order(abs(pvalue), decreasing=FALSE)
lpv = ordered_lpv[1:length(pvalue[np])] 

oo = union(lpv, hfc) 
ii = intersect(lpv, hfc) 


## @knitr F3
#Construct a volcano plot using moderated F-statistics
plot(difference[-oo], lodd[-oo], xlim=range(difference), 
ylim=range(lodd), cex=0.7)
points(difference[hfc], lodd[hfc], pch=18, cex=0.7) 
points(difference[lpv], lodd[lpv], pch=1, cex=0.7)

#Add the cut-off lines
abline(v=log(2, 2), col=5); abline(v=-log(2, 2), col=5)
abline(h=-log10(0.05/5900), col=5)

text(min(difference) + 1, -log10(0.05/5900) + 0.2, 
    'Bonferroni cut off', cex=0.8) 
text(1, max(lodd) - 1, paste(length(ii), 'intersects'), cex=0.8)


## @knitr unnamed-chunk-23
c_probe_data = yeast.matrix[ii,]
#Average of WT
wt_means = apply(c_probe_data[,16:30], 1, mean)
m = matrix(nrow=dim(c_probe_data)[1], ncol=5)

for(i in 1:5){
    mut_rep = c(i, i+5, i+10)
    m[ ,i] = rowMeans(c_probe_data[ ,mut_rep]) - wt_means
}
colnames(m) = sort(unique(exp_fac$tps))


## @knitr F4
library(gplots)
#Cluster the top 50 genes
heatmap.2(m[1:50,], dendrogram ='row', Colv=FALSE, col=greenred(75), 
          key=FALSE, keysize=1.0, symkey=FALSE, density.info='none', 
          trace='none', colsep=rep(1:10), sepcolor='white', 
          sepwidth=0.05, labRow = NA, cexCol=1,
          hclustfun=function(c){hclust(c, method='average')})


## @knitr F5
library(Mfuzz)
tmp_expr = new('ExpressionSet', exprs=m)
cl = mfuzz(tmp_expr, c=8, m=1.25)
mfuzz.plot(tmp_expr,cl=cl, mfrow=c(2, 4), new.window = FALSE)


## @knitr unnamed-chunk-24
cluster = 1
cl[[4]][, cluster] 


## @knitr unnamed-chunk-25
exp_fac = with(exp_fac, exp_fac[order(strain, tps, replicates), ])
#Construct a longitudinal object
library(GeneNet)
ngenes = 100
m = yeast.matrix[ii[1:ngenes],]
mnew = m[,exp_fac$data_order[1:15]]
mlong = as.longitudinal(t(mnew), repeats=3, time=0:4)


## @knitr unnamed-chunk-26
#Compute partial correlations
pcor.dyn = ggm.estimate.pcor(mlong, method = 'dynamic')
 
#Assign (local) fdr values to all possible edges
m.edges = network.test.edges(pcor.dyn, direct=TRUE)

#Construct graph containing top edges
m.net = extract.network(m.edges, method.ggm='number', 
cutoff.ggm=100)

#Construct a Graphviz dot file
rnames = vector('list', length(1))
rnames = c_df$genename[ii[1:ngenes]]
network.make.dot(filename='net.dot', m.net, rnames, 
    main='Yeast Network')


## @knitr ExtractIDs



## @knitr RemoveProbes



## @knitr unnamed-chunk-27
#Function to average the expression of 
#probesets which map to same gene
probeset2genelevel = function(onesample)
    return(tapply(onesample, factor(c_df$genename), mean))

#Average for each column/array
c_gene_data = apply(exprs(yeast.rma), 2, probeset2genelevel)


## @knitr F1A
par(mar=c(3,3,2,1), 
           mgp=c(2,0.4,0), tck=-.01,
               cex.axis=0.9, las=1, mfrow=c(2,3))
run_chunk("s_fig1")


## @knitr F2A
par(mar=c(3,3,2,1), 
           mgp=c(2,0.4,0), tck=-.01,
               cex.axis=0.9, las=1)
run_chunk("s_fig2")


## @knitr F3A
par(mar=c(3,3,2,1), 
           mgp=c(2,0.4,0), tck=-.01,
               cex.axis=0.9, las=1, mfrow=c(1,2))
run_chunk("s_fig3")


## @knitr F4A
par(mar=c(3,3,2,1), 
           mgp=c(2,0.4,0), tck=-.01,
               cex.axis=0.9, las=1)
run_chunk("s_fig4")


## @knitr unnamed-chunk-28
sessionInfo()


