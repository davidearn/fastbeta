SHELL = /bin/bash
name = fastbeta-vignette

all : $(name).pdf

$(name).pdf : $(name).tex $(name).bib
	pdflatex -interaction=batchmode $<
	bibtex $(name)
	pdflatex -interaction=batchmode $<
	pdflatex $<

$(name).tex: $(name).Rnw knit_theme.css
	Rscript -e "knitr::knit('$<')"

clean :
	rm -rf cache figure
	rm -f *.{bbl,log,synctex.gz,toc,aux,tex}
