bmc-microarray
==============

This repository contains the source code for the paper: Gillespie, C. S., et al, 2010. [Analysing yeast time course microarray data using BioConductor: a case study using yeast2 Affymetrix arrays](http://www.biomedcentral.com/1756-0500/3/81). *BMC Research Notes*, **3:81**. Since this paper is now a few years old, the original code does not work. This repository contains a slightly modified version of the paper to work with the latest versions of R and Bioconductor. 

### Documents

 * The updated paper ([pdf](https://github.com/csgillespie/bmc-microarray/raw/master/paper.pdf))
 * [R code](https://github.com/csgillespie/bmc-microarray/raw/master/paper.R) from the paper
 * Additional functions: [ExtractIDs.R](https://github.com/csgillespie/bmc-microarray/raw/master/ExtractIDs.R), [RemoveProbes.R](https://github.com/csgillespie/bmc-microarray/raw/master/RemoveProbes.R)
 * Mask file: s_cerevisiae.msk ([direct link](https://github.com/csgillespie/bmc-microarray/raw/master/s_cerevisiae.msk), [affymetrix web site](http://www.affymetrix.com/Auth/support/downloads/mask_files/s_cerevisiae.zip))
 * Annotation file: yeast2annotation.csv ([direct link](https://github.com/csgillespie/bmc-microarray/raw/master/yeast2annotation.csv), [affymetrix web site](http://www.affymetrix.com/Auth/analysis/downloads/na24/ivt/Yeast_2.na24.annot.csv.zip))
 

#### Updates/Corrections since publication

  * (11/07/2012) Changed from Sweave to knitr. Moved source to github.
  * (10/07/2012) The *yeast2probe* environment became protected. Minor changes to the *RemoveProbes.R* function. Thanks to Guang You Duan for pointing this out.
  * (11/06/2010) Paper updated to reflect the new array express download format. Thanks to Saeed Salem for pointing this out.

#### Suggestions, bugs, queries

If you have any questions or queries about the above material, then feel free to contact me at csgillespie@gmail.com Also if you think of any other techniques that could be included in the paper, then please send in your suggestions.

#### Building the paper from source

The paper was originally constructed in Sweave (thanks to a referees comment!). It has now been ported over to knitr. To build the pdf of the paper, just run the following commands (this assumes you have the necessary cran and Bioconductor libraries and have downloaded all files from github):

```{r}
library(knitr)
fname = 'paper.Rnw'
knit(fname)
## extract R code only
purl(fname)
```

Alternatively, you can use an IDE such as RStudio.

