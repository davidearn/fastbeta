SHELL = bash
R = R
PACKAGE = fastbeta
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
VIGNETTE = fastbeta-vignette.Rnw fastbeta-vignette.bib knit_theme.css Makefile

basic :
	$(R) --no-echo -e 'devtools::document(".")'
	$(R) CMD build --no-manual .
	$(R) CMD INSTALL $(TARBALL)

fancy : rd namespace install

all : clean install check

install : build
	export NOT_CRAN=true
	$(R) CMD INSTALL ../$<

build : $(TARBALL)

$(TARBALL) : R/*.R man/*.Rd vignettes/$(VIGNETTE) DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .
	mv $@ ..

man/*.Rd &: R/*.R inst/REFERENCES.bib
	$(R) --no-echo -e 'devtools::document(".", roclets = "rd")'
	@touch $@

DESCRIPTION : R/*.R
	$(R) --no-echo -e 'devtools::document(".", roclets = "collate")'
	@touch $@

NAMESPACE : R/*.R
	$(R) --no-echo -e 'devtools::document(".", roclets = "namespace")'
	@touch $@

check :
	$(R) --no-echo -e 'devtools::check(".")'

test :
	$(R) --no-echo -e 'devtools::test(".")'

clean :
	rm -f ../$(TARBALL)
	find . \( -name "\.#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +
