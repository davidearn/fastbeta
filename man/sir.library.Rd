\name{sir.library}
\alias{sir.library}
\alias{sir.e01}
\alias{sir.e02}
\docType{data}
\title{Often Used Simulations}
\description{
Infectious disease time series simulated using \code{\link{sir}},
for use primarily in examples, tests, and vignettes.  Users should
not rely on simulation details, which may change between package
versions.
}
\usage{
data(sir.e01, package = "fastbeta") # "Example 01", and so on
data(sir.e02, package = "fastbeta")
}
\format{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, always a subset of the result of a call
to \code{\link{sir}}, discarding transient behaviour.
Simulation parameters may be preserved as attributes.
}
\source{
Scripts sourced by \code{\link{data}} to reproduce the simulations
are located in the \file{data} subdirectory of the \pkg{fastbeta}
installation; see, e.g.
\code{\link{system.file}("data", "sir.e01.R", package = "fastbeta")}.
}
\seealso{\code{\link{sir}}.}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
data(sir.e01, package = "fastbeta")
str(sir.e02)
plot(sir.e01)

data(sir.e02, package = "fastbeta")
str(sir.e02)
plot(sir.e02)
}
