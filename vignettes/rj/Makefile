.PHONY : all
all : RJwrapper.pdf

RJwrapper.pdf : RJwrapper.tex fastbeta-ms.tex fastbeta-ms.bib
	R --quiet --vanilla -e 'tools::texi2pdf("$<")'

fastbeta-ms.tex : fastbeta-ms.Rnw
	R --quiet --vanilla -e 'utils::Sweave("$<")'

.PHONY : clean
clean :
	rm -f fastbeta-ms.tex
	rm -f *~ .*~ *.aux *.bbl *.brf *.blg *.log *.out

.PHONY : distclean
distclean : clean
	rm -f RJwrapper.pdf
