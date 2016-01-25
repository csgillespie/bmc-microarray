# Makefile

## $* = filename without extension
## $@ = the output file
## $< = the input file
.SUFFIXES: .tex .pdf .Rnw .R
PAPER=paper

.PHONY: clean force

$(PAPER).pdf: $(PAPER).tex references.bib

.Rnw.tex:
	Rscript -e 'fname = "paper.Rnw"; knitr::knit(fname); knitr::purl(fname)'

.tex.pdf:
	echo pwd
	pdflatex $*.tex
	bibtex $*
	pdflatex $*.tex
	pdflatex $*.tex


pdf:
	make $(PAPER).pdf


R: *.Rnw
	make cleanR	
	R CMD Stangle $(PAPER)
	mv ExtractIDs*.R ExtractIDs.R
	mv RemoveProbes*.R RemoveProbes.R
#	#Remove columns and first few rows
	colrm 1 2 < $(PAPER).R > tmp.R
	tail -n+2 tmp.R > $(PAPER).R

#	#TODO: Learn to use Makefiles correctly
	colrm 1 2 < ExtractIDs.R > tmp.R; mv tmp.R ExtractIDs.R
	colrm 1 2 < RemoveProbes.R > tmp.R; mv tmp.R RemoveProbes.R



all: *.Rnw refs.bib 
	make R
	make $(PAPER).pdf


clean:
	rm -f $(PAPER).ps $(PAPER).pdf $(PAPER).tex $(PAPER).R  tmp.R
	rm -f *~ *.aux *.toc *.out *.log *.dvi *.bak *.snm *.nav *.vrb *.bbl *.blg 
	rm -fv updown.csv genelist.csv indexE-MEXP-1551.html
	rm -f  net.dot

# eof

