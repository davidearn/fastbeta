package := $(shell sed -n 's/Package: //p' DESCRIPTION)
version := $(shell sed -n 's/Version: //p' DESCRIPTION)
tarball := $(package)_$(version).tar.gz
rchkdir := $(package)_$(version).Rcheck

sources := .Rbuildignore DESCRIPTION NAMESPACE \
	R/*.R man/*.Rd tests/*.R vignettes/*.Rnw vignettes/*.bib vignettes/*.tex

all: build

build: $(tarball)

$(tarball): $(sources)
	R CMD build --compact-vignettes=gs+qpdf .

install: $(tarball)
	R CMD INSTALL --preclean --clean $<

check: $(tarball)
	R CMD check -o $(rchkdir) $<

clean:
	rm -rf $(tarball) $(rchkdir)
	find . -name '*~' -type f -exec rm {} \+
