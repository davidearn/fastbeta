SHELL = bash
R = R
PACKAGE = fastbeta
VERSION := $(shell sed -n '/^Version: /s///p' ./DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz

basic:
	Rscript -e "devtools::document()"
	R CMD build --no-manual .
	R CMD INSTALL *.tar.gz

fancy: doc-update install

pkgall: clean doc-update namespace-update install pkgcheck

doc-update: R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\".\",roclets = c(\"collate\", \"rd\")))" | $(R) --slave
	@touch $@

namespace-update: R/*.R
	echo "suppressWarnings(roxygen2::roxygenize('.',roclets = 'namespace'))" | $(R) --slave
	@touch $@

$(TARBALL): ./NAMESPACE
	$(R) CMD build .
	mv $@ ..

install: $(TARBALL)
	##export NOT_CRAN=true; $(R) CMD INSTALL --preclean --no-manual ../$<
	$(R) CMD INSTALL *.tar.gz
	@touch $@

pkgtest:
	echo "devtools::test('.')" | $(R) --slave

pkgcheck:
	echo "devtools::check('.')" | $(R) --slave

clean:
	rm -f *.tar.gz
	find . \( -name "\.#*" -o -name "*~" -o -name ".Rhistory" \) -exec rm {} \;
