bmc-microarray
==============

This repository contains the source code for the paper: Gillespie, C. S., et al, 2010. [Analysing yeast time course microarray data using BioConductor: a case study using yeast2 Affymetrix arrays](http://www.biomedcentral.com/1756-0500/3/81). *BMC Research Notes*, **3:81**. Since this paper is now a few years old, the original code does not work. This repository contains a slightly modified version of the paper to work with the latest versions of R and Bioconductor. 

The paper was originally constructed in Sweave (thanks to a referees comment!). It has now been ported over to knitr. To build the pdf of the paper, just run the following commands (this assumes you have the necessary cran and Bioconductor libraries):

```{r}
library(knitr)
fname = 'paper.Rnw'
knit(f)
## extract R code only
purl(f)
```

Alternatively, you can use an IDE such as RStudio.

Occasionally, the code in the paper breaks due to an update in R or bioconductor. If this happens, please contact me.
