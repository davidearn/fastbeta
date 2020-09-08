SHELL = /bin/bash
R_HOME = /usr/lib/R
name = fastbeta-vignette

all : $(name).pdf new

$(name).pdf : $(name).tex $(name).bib
	pdflatex -interaction=batchmode $<
	bibtex $(name)
	pdflatex -interaction=batchmode $<
	pdflatex $<

$(name).tex: $(name).Rnw knit_theme.css
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')"

new : $(name).pdf
	rm -rf cache figure
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc,aux}

clean :
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc}
