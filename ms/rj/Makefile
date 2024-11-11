R := R --quiet --vanilla
ms := fastbeta-ms
main := RJwrapper

.PHONY : all
all : $(ms).pdf

$(ms).pdf : $(main).tex $(ms).tex $(ms).bib
	$(R) -e 'tools::texi2pdf("$(main).tex")'
	cp $(main).pdf $(ms).pdf

$(ms).tex : $(ms).Rnw
	$(R) -e 'utils::Sweave("$(ms).Rnw")'
	sed \
	    -i '' \
	    -e 's/\begin{Schunk}/\begin{example}/g' \
	    -e 's/\end{Schunk}/\end{example}/g' \
	    -e '/\begin{Sinput}/d' \
	    -e '/\end{Sinput}/d' \
	    -e '/\begin{Soutput}/d' \
	    -e '/\end{Soutput}/d' \
	    $@

.PHONY : clean
clean :
	rm -f $(ms).tex figure-[0-9]*.tex
	rm -f *~ .*~ *.aux *.bbl *.brf *.blg *.log *.out

.PHONY : distclean
distclean : clean
	rm -f $(ms).pdf $(main).pdf figure-[0-9]*.pdf
