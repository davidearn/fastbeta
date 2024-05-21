# fastbeta

**fastbeta** is an R package for approximating time-varying infectious
disease transmission rates from disease incidence time series and other
data.  The algorithm is fast, being based on a discrete time
approximation of an SEIR model with user-defined numbers of latent and
infectious compartments.  Auxiliary functions support parametric
bootstrapping (for approximating confidence intervals on transmission
rates), efficient time series simulation (for general purpose simulation
studies), and Richardson-Lucy deconvolution (for estimating incidence
from *observed* incidence conditional on a flexible observation model).

## Installation

[CRAN](https://cran.r-project.org/package=fastbeta) distributes
both the package sources and binaries for Windows and macOS.  Hence
typical users will install **fastbeta** with

```r
install.packages("fastbeta")
```

or perhaps

```r
install.packages("fastbeta", type = "source")
```

to force installation from sources where installation of a binary
would occur by default.  The rest of this section concerns the
`type = "source"` case.

Installation from sources depends on compilers and related tools.
These will already be available on modern Linux installations.
Windows users must have installed
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).
macOS users must have installed Apple's Command Line Tools for
[Xcode](https://developer.apple.com/xcode/)
and GNU Fortran.
The most recent version of Command Line Tools supporting the
version of macOS in use can be installed by running

```shell
sudo rm -rf /Library/Developer/CommandLineTools
sudo xcode-select --install
```

in Terminal.
Binaries for older versions of Command Line Tools can be downloaded
[here](https://developer.apple.com/download/all/?q=Command%20Line%20Tools%20for%20Xcode).
GNU Fortran should be installed following the instructions on the
R for macOS Developers [page](https://mac.r-project.org/tools/).

Vignette building depends on a
[LaTeX](https://www.latex-project.org/get/) distribution.
Specifically, `PATH` must specify the location of `pdflatex`.

Issues related to installation should be reported (with relevant output)
[here](https://github.com/davidearn/fastbeta/issues/1).
Notably, details for Windows and macOS can differ for non-current
versions of R or non-standard installations of R, where
the "standard" way to install R is to download and unpack a binary
built and published by [CRAN](https://cran.r-project.org/).

Compiler errors encountered on macOS are almost always explained
by unmet dependencies or masking of native tools, headers, and
libraries with non-native ones (e.g., ones installed by Homebrew).
Masking occurs due to dubious configuration of `PATH` or dubious
setting (typically in `~/.R/Makevars`) of Make variables such as
`CPPFLAGS` and `LDFLAGS`.  Users should reattempt compilation
after removing suspicious components of `PATH` (e.g., by
removing relevant lines of startup files in your home directory,
then launching a new shell) and (re)moving `~/.R/Makevars`.

## Documentation

After installing, users can access the package index (a list of
available help topics) with:

```r
help(package = "fastbeta")
```

The HTML help contains useful hyperlinks and typeset math.
You can force HTML rendering, where that is not the default,
by passing `help_type = "html"`.

## Repository structure

Active development happens on branch `master`.  Tested changes intended
for the next release are ported to branch `release-candidate`,
where tarballs submitted to CRAN are eventually built.  Neither `master`
nor `release-candidate` should be considered stable.

The stable branches are named `release-x.y.z`.  They branch from
`release-candidate` before the version number there is incremented,
typically just after a tarball is submitted to CRAN.

To install **fastbeta** from sources in a given branch or commit,
install [**remotes**](https://cran.r-project.org/package=remotes) and
run, e.g.,

```r
remotes::install_github("davidearn/fastbeta", ref = "release-0.3.0")
remotes::install_github("davidearn/fastbeta", ref = "04c8f7b")
```

## References

Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
Fast estimation of time-varying infectious disease transmission rates.
*PLOS Computational Biology*, *16*(9), Article e1008124, 1-39.
[https://doi.org/10.1371/journal.pcbi.1008124]
