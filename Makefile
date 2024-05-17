package := $(shell sed -n 's/Package: //p' DESCRIPTION)
version := $(shell sed -n 's/Version: //p' DESCRIPTION)
tarball := $(package)_$(version).tar.gz
rchkdir := $(package).Rcheck
rchklog := $(rchkdir)/00check.log

sources := \
	.Rbuildignore \
	DESCRIPTION \
	NAMESPACE \
	R/*.R \
	data/*.R \
	data/*.rda \
	data/datalist \
	inst/NEWS.Rd \
	inst/scripts/*.R \
	man/*.Rd \
	src/*.c \
	tests/*.R

.PHONY: all
all: build

.PHONY: build
build: $(tarball)

.PHONY: check
check: $(rchklog)

.PHONY: install
install: $(tarball)
	R CMD INSTALL $<

.PHONY: clean
clean:
	rm -rf $(tarball) $(rchkdir)
	find . -name *~ -type f -exec rm {} \+

$(tarball): $(sources)
	R CMD build --no-build-vignettes .

$(rchklog): $(tarball)
	R CMD check --as-cran $<
