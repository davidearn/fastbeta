.PHONY : all
all : RJwrapper.pdf

RJwrapper.pdf : RJwrapper.tex fastbeta-ms.tex fastbeta-ms.bib
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

fastbeta-ms.tex : fastbeta-ms.Rnw
	R --vanilla -e 'Sweave("$<")'

.PHONY : clean
clean :
	rm -f fastbeta-ms.tex
	rm -f *~ .*~ *.aux *.bbl *.brf *.blg *.log *.out

.PHONY : distclean
distclean : clean
	rm -f RJwrapper.pdf
