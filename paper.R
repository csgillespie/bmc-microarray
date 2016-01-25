## ----echo=FALSE----------------------------------------------------------
library("knitr")
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

## ----eval=FALSE----------------------------------------------------------
## url = 'http://bioconductor.org/biocLite.R'
## source(url)
## biocLite()

## ----eval=FALSE,tidy=FALSE-----------------------------------------------
## ##From Bioconductor
## biocLite(c('ArrayExpress', 'Mfuzz', 'timecourse', 'yeast2.db',
##     'yeast2probe', 'yeast2cdf', 'affyPLM'))
## ##From cran
## install.packages(c('GeneNet', 'gplots'))

## ----eval=FALSE----------------------------------------------------------
## update.packages(repos=biocinstallRepos())

## ----eval=TRUE,hide=TRUE,cache=TRUE,message=FALSE------------------------
library('ArrayExpress')
yeast.raw = ArrayExpress('E-MEXP-1551')

## ----eval=TRUE,cache=TRUE,tidy=FALSE-------------------------------------
ph = yeast.raw@phenoData
exp_fac = data.frame(data_order = 1:30, 
    strain = ph@data$Factor.Value..phenotype.,  
    replicates = ph@data$Characteristics..individual.,
    tps = ph@data$Factor.Value..time.)
levels(exp_fac$strain) = c('m', 'w')
exp_fac = with(exp_fac, exp_fac[order(strain, replicates, tps), ])
exp_fac$replicate = rep(1:3, each=5, times=2)

## ----eval=TRUE,message=FALSE,warning=FALSE-------------------------------
#Read in the mask file
s_cer = read.table('s_cerevisiae.msk', skip=2, 
    stringsAsFactors=FALSE)
source('ExtractIDs.R')
c_df = ExtractIDs(s_cer[ ,1])

## ----eval=TRUE,message=FALSE---------------------------------------------
#Get the raw dataset for S. cerevisiae only
library('affy')
library('yeast2probe')
source('RemoveProbes.R')
cleancdf = cleancdfname(yeast.raw@cdfName)
RemoveProbes(s_cer[, 1], cleancdf, 'yeast2probe') 

## ----s_fig1,eval=FALSE,tidy=TRUE, echo=-1--------------------------------
## par(mar=c(3,3,2,1), mgp=c(2,0.4,0), tck=-.01, cex.axis=0.9, las=1, mfrow=c(2,3))
## for(i in 1:5) {
##     plot_title = paste('Strain: ', exp_fac$strain[i],
##         'Time: ', exp_fac$tps[i])
##         d = exp_fac$data_order[i]
##     image(yeast.raw[,d], main=plot_title)
## }

## ----s_fig2,eval=FALSE,tidy=FALSE,echo=-1--------------------------------
## par(mar=c(3,3,2,1), mgp=c(2,0.4,0), tck=-.01, cex.axis=0.9, las=1)
## d = exp_fac$data_order[1:5]
## hist(yeast.raw[,d], lwd=2, ylab="Density", xlab="Log (base 2) intensities")

## ----eval=TRUE,cache=TRUE,results=FALSE,message=FALSE--------------------
yeast.rma = rma(yeast.raw)
yeast.matrix = exprs(yeast.rma)[, exp_fac$data_order]
colnames(yeast.matrix) = paste0(exp_fac$strain, exp_fac$tps)
exp_fac$data_order = 1:30

## ----s_fig3,eval=FALSE, echo=-1------------------------------------------
## par(mar=c(3,3,2,1), mgp=c(2,0.4,0), tck=-.01,cex.axis=0.9, las=1, mfrow=c(1,2))
## library('affyPLM')
## #Raw data intensities
## boxplot(yeast.raw, col='red', main='', ylim=c(2, 16))
## #Normalised intensities
## boxplot(yeast.rma, col='blue', ylim=c(2, 16))

## ----eval=TRUE-----------------------------------------------------------
yeast.PC = prcomp(t(yeast.matrix))
yeast.scores = predict(yeast.PC)

## ----F1,eval=TRUE, tidy=FALSE,dev="pdf",fig.pos="!t", cache=TRUE, fig.width=6.5,fig.height=4,out.width="0.8\\textwidth",fig.cap="A plot of the first two principal components. The red symbols correspond to the wild-type strain.", par.nice=TRUE----
#Plot of the first two principal components
plot(yeast.scores[,1], yeast.scores[,2], xlab='PC 1', ylab='PC 2', 
    pch=rep(1:5, 6), col=as.numeric(exp_fac$strain))
legend("bottomleft", pch=1:5, cex=0.6, 
    c('t 0', 't 60', 't 120', 't 180', 't 240'))

## ----eval=TRUE,message=FALSE---------------------------------------------
library('timecourse')
size = matrix(3, nrow = 5900, ncol = 2)

## ----eval=TRUE,tidy=FALSE, cache=TRUE------------------------------------
c.grp = as.character(exp_fac$strain)
t.grp = as.numeric(exp_fac$tps)
r.grp = as.character(exp_fac$replicate)
MB.2D = mb.long(yeast.matrix, times = 5, method = '2', 
    reps = size, condition.grp = c.grp, time.grp = t.grp, 
    rep.grp = r.grp)

## ----eval=TRUE,cache=TRUE------------------------------------------------
gene_positions = MB.2D$pos.HotellingT2[1:100]
gnames = rownames(yeast.matrix)
gene_probes = gnames[gene_positions]

## ----s_fig4,eval=FALSE, echo=-1------------------------------------------
## par(mar=c(3,3,2,1), mgp=c(2,0.4,0), tck=-.01, cex.axis=0.9, las=1)
## plotProfile(MB.2D, ranking=1, gnames=rownames(yeast.matrix))

## ----eval=TRUE,tidy=FALSE, message=FALSE---------------------------------
library('limma')
expt_structure = factor(colnames(yeast.matrix))

#Construct the design matrix
X = model.matrix(~0 + expt_structure)
colnames(X) =  c('m0', 'm60', 'm120', 'm180', 'm240', 
    'w0', 'w60', 'w120', 'w180', 'w240')

## ----eval=TRUE,cache=TRUE------------------------------------------------
lm.fit = lmFit(yeast.matrix, X) 

## ----eval=TRUE,cache=TRUE------------------------------------------------
mc = makeContrasts('m60-w60', 'm120-w120', 'm180-w180', 'm240-w240', levels=X)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)

## ----eval=TRUE,results='hide',cache=TRUE---------------------------------
#see help(toptable) for more options
toptable(eb, sort.by='logFC') 

## ----eval=TRUE, results='hide',cache=TRUE--------------------------------
topTableF(eb)

## ----eval=TRUE,cache=TRUE,tidy=FALSE-------------------------------------
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

## ----F2,eval=TRUE,echo=FALSE,dev="pdf",fig.pos="!t", cache=TRUE, fig.width=6.5,fig.height=5,out.width="0.8\\textwidth",fig.cap="Time course expression levels for the top 9 differentially expressed genes, ranked by their $F$-statistic. The triangles and circles correspond to the wild-type and mutant genes respectively.",par.nice=TRUE, par.mfrow=c(3,3), tidy=FALSE----

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

## ----eval=FALSE, tidy=FALSE----------------------------------------------
## #Rank of Probesets, also output gene names
## par(mfrow=c(3, 3), ask=TRUE)
## for(i in 0:99){
##     indx = rank(modF) == nrow(yeast.matrix) - i
##     id = c_df$probe[indx]
##     name = c_df$genename[indx]
##     exprs.row = yeast.matrix[indx,]
##     genetitle = paste(sprintf('%.30s', id), sprintf('%.30s', name),
##                       'Rank =', i+1)
## 
##     plot(0, pch=NA, xlim=range(0, 240), ylim=range(exprs.row),
##          ylab='Expression', xlab='Time', main=genetitle)
## 
##     for(j in 1:6){
##         pch_value = as.character(exp_fac$strain[5*j])
##         points(c(0, 60, 120, 180, 240), exprs.row[(5*j-4):(5*j)],
##                type='b', pch=pch_value)
##     }
## }

## ----eval=TRUE, cache=TRUE-----------------------------------------------
N = 100
gene_positions = MB.2D$pos.HotellingT2[1:N]
tc_top_probes = gnames[gene_positions]
lm_top_probes = c_df$probe[modFordered[1:N]]
length(intersect(tc_top_probes, lm_top_probes))

## ----eval=TRUE,cache=TRUE------------------------------------------------
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

## ----F3,eval=TRUE,dev="pdf",fig.pos="!t", cache=TRUE, fig.width=6.5,fig.height=4,out.width="0.8\\textwidth",fig.cap="Volcano plot showing the bonferroni cut-off and the two-fold change.", nice.par=TRUE,tidy=FALSE----
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

## ----eval=TRUE-----------------------------------------------------------
c_probe_data = yeast.matrix[ii,]
#Average of WT
wt_means = apply(c_probe_data[,16:30], 1, mean)
m = matrix(nrow=dim(c_probe_data)[1], ncol=5)

for(i in 1:5){
    mut_rep = c(i, i+5, i+10)
    m[ ,i] = rowMeans(c_probe_data[ ,mut_rep]) - wt_means
}
colnames(m) = sort(unique(exp_fac$tps))

## ----F4,eval=TRUE,dev="pdf",fig.pos="!t", cache=TRUE, fig.width=6.5,fig.height=4,out.width="0.8\\textwidth",fig.cap="Clustering of the top fifty differentially expressed genes. Red and green correspond to up- and down-regulation respectively.",message=FALSE, nice.par=TRUE,tidy=FALSE----
library('gplots')
# Cluster the top 50 genes
heatmap.2(m[1:50,], dendrogram ='row', Colv=FALSE, col=greenred(75), 
          key=FALSE, keysize=1.0, symkey=FALSE, density.info='none', 
          trace='none', colsep=rep(1:10), sepcolor='white', 
          sepwidth=0.05, labRow = NA, cexCol=1,
          hclustfun=function(c){hclust(c, method='average')})

## ----F5,eval=TRUE,dev="pdf",fig.ext='pdf',fig.pos="!t", cache=TRUE, fig.width=6.5,fig.height=4,out.width="0.8\\textwidth",fig.cap="Eight clusters obtained from the \\texttt{Mfuzz} package.",message=FALSE----
library('Mfuzz')
tmp_expr = new('ExpressionSet', exprs=m)
cl = mfuzz(tmp_expr, c=8, m=1.25)
mfuzz.plot(tmp_expr,cl=cl, mfrow=c(2, 4), new.window = FALSE)

## ----eval=TRUE,results='hide'--------------------------------------------
cluster = 1
cl[[4]][, cluster] 

## ----eval=TRUE,results='hide',message=FALSE------------------------------
exp_fac = with(exp_fac, exp_fac[order(strain, tps, replicates), ])
#Construct a longitudinal object
library('GeneNet')
ngenes = 100
m = yeast.matrix[ii[1:ngenes],]
mnew = m[,exp_fac$data_order[1:15]]
mlong = as.longitudinal(t(mnew), repeats=3, time=0:4)

## ----eval=FALSE,cache=TRUE,message=FALSE,results='hide'------------------
## #Compute partial correlations
## pcor.dyn = ggm.estimate.pcor(mlong, method = 'dynamic')
## 
## #Assign (local) fdr values to all possible edges
## m.edges = network.test.edges(pcor.dyn, direct=TRUE)
## 
## #Construct graph containing top edges
## m.net = extract.network(m.edges, method.ggm='number',
## cutoff.ggm=100)
## 
## #Construct a Graphviz dot file
## rnames = vector('list', length(1))
## rnames = c_df$genename[ii[1:ngenes]]
## network.make.dot(filename='net.dot', m.net, rnames,
##     main='Yeast Network')

## ----ExtractIDs,tidy=FALSE-----------------------------------------------

## ----RemoveProbes,tidy=FALSE---------------------------------------------

## ----eval=FALSE----------------------------------------------------------
## #Function to average the expression of
## #probesets which map to same gene
## probeset2genelevel = function(onesample)
##     return(tapply(onesample, factor(c_df$genename), mean))
## 
## #Average for each column/array
## c_gene_data = apply(exprs(yeast.rma), 2, probeset2genelevel)

## ----F1A,ref.label="s_fig1",eval=TRUE,dev="png",fig.ext="png",cache=FALSE,dpi=300,out.width="0.8\\linewidth",echo=FALSE, fig.cap="Image plots of the mismatch and perfect match probe intensities for the first replication of the mutant yeast strain. The corresponding times are indicated in the plot.", fig.pos="h"----
#knitr::run_chunk("s_fig1")

## ----F2A,ref.label="s_fig2",eval=TRUE,dev="pdf",fig.ext="pdf",cache=TRUE,out.width="0.6\\linewidth",echo=FALSE, fig.cap="Density plots for the first replication of the mutant yeast strain.",nice.par=TRUE,message=FALSE----

## ----F3A,ref.label="s_fig3",eval=TRUE,dev="pdf",fig.ext="pdf",cache=TRUE,out.width="0.8\\linewidth",echo=FALSE, message=FALSE,fig.cap="Boxplots of the raw and normalised intensities. The default boxplot is to include both PM and MM intensities, whereas for the density plots in Figure~\ref{F2} the default is for only the PM intensities."----

## ----F4A,ref.label="s_fig4", eval=TRUE,dev="pdf",fig.ext="pdf",cache=TRUE,out.width="0.6\\linewidth",echo=FALSE, fig.cap="Time course expression levels for the top differentially expressed gene, ranked by their Hotelling statistic using the \texttt{timecourse} library."----

## ----eval=TRUE-----------------------------------------------------------
sessionInfo()

