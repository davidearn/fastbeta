SHELL = /bin/bash

all :
	Rscript -e "devtools::document()"
	R CMD build --no-manual .
	R CMD INSTALL *.tar.gz
	rm -f *.tar.gz
