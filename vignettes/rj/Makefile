R := R --quiet --vanilla

  name := fastbeta-ms
rjname := RJwrapper

.PHONY : all
all : $(rjname).pdf

$(rjname).pdf : $(rjname).tex $(name).tex $(name).bib
	$(R) -e 'tools::texi2pdf("$<")'

$(name).tex : $(name).Rnw
	$(R) -e 'utils::Sweave("$<")'
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
	rm -f $(name).tex
	rm -f *~ .*~ *.aux *.bbl *.brf *.blg *.log *.out

.PHONY : distclean
distclean : clean
	rm -f $(rjname).pdf
