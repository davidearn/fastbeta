package := $(shell sed -n 's/Package: //p' DESCRIPTION)
version := $(shell sed -n 's/Version: //p' DESCRIPTION)
tarball := $(package)_$(version).tar.gz
rchkdir := $(package).Rcheck
rchklog := $(rchkdir)/00check.log

sources := .Rbuildignore DESCRIPTION NAMESPACE \
	R/*.R data/*.rda man/*.Rd tests/*.R \
	vignettes/*.Rnw vignettes/*.bib vignettes/*.tex

.PHONY: all build check install clean

all: build

build: $(tarball)

check: $(rchklog)

install: $(tarball)
	R CMD INSTALL $<

clean:
	rm -rf $(tarball) $(rchkdir)
	find . -name '*~' -type f -exec rm {} \+

$(tarball): $(sources)
	R CMD build --no-resave-data --compact-vignettes=gs+qpdf .

$(rchklog): $(tarball)
	R CMD check --as-cran $<
