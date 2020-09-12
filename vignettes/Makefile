SHELL = bash
R = R
NAME = fastbeta-vignette

build: $(NAME).pdf
	rm -rf cache figure
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc,aux}

$(name).pdf: $(NAME).tex $(NAME).bib
	pdflatex -interaction=batchmode $<
	bibtex $(NAME)
	pdflatex -interaction=batchmode $<
	pdflatex $<

$(name).tex: $(NAME).Rnw knit_theme.css
	$(R) --no-echo -e 'knitr::knit("$<")'

all: $(name).pdf

new: clean
	rm -rf cache figure *.{aux,pdf}

clean:
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc}
