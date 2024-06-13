.PHONY : all
all : RJwrapper.pdf

RJwrapper.pdf : RJwrapper.tex fastbeta-ms.tex fastbeta-ms.bib
	R --quiet --vanilla -e 'tools::texi2pdf("$<")'

fastbeta-ms.tex : fastbeta-ms.Rnw
	R --quiet --vanilla -e 'utils::Sweave("$<")'
	sed -i '' -e 's/\begin{Schunk}/\begin{example}/g' \
	          -e 's/\end{Schunk}/\end{example}/g' \
	          -e '/\begin{Sinput}/d' \
	          -e '/\end{Sinput}/d' \
	          -e '/\begin{Soutput}/d' \
	          -e '/\end{Soutput}/d' $@

.PHONY : clean
clean :
	rm -f fastbeta-ms.tex
	rm -f *~ .*~ *.aux *.bbl *.brf *.blg *.log *.out

.PHONY : distclean
distclean : clean
	rm -f RJwrapper.pdf
