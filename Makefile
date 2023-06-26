package := $(shell grep DESCRIPTION -e Package: | sed 's/Package: //')
version := $(shell grep DESCRIPTION -e Version: | sed 's/Version: //')
tarball := $(package)_$(version).tar.gz
rchkdir := $(package)_$(version).Rcheck

all: build

build: $(tarball)

$(tarball):
	R CMD build --compact-vignettes=gs+qpdf .

install: $(tarball)
	R CMD INSTALL --preclean --clean $<

check: $(tarball)
	R CMD check -o $(rchkdir) $<

clean:
	rm -rf $(tarball) $(rchkdir)
	find . -name '*~' -type f -exec rm {} \+
