\name{sir.library}
\alias{sir.library}
\alias{sir.e01}
\docType{data}
\title{Often used simulations}
\description{
Infectious disease time series simulated using \code{\link{sir}},
for use primarily in examples, tests, and vignettes.
}
\usage{
   data(sir.e01, package = "fastbeta") # "Example 01"
## data(sir.e02, package = "fastbeta") # "Example 02" (TODO)
}
\format{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, always a subset of the result of a call
to \code{\link{sir}}, discarding transient behaviour.
Simulation parameters may be preserved as attributes.
}
\source{
To reproduce these objects, run the corresponding scripts, e.g.,
\code{\link{system.file}("data", "sir.e01.R", package = "fastbeta")}.
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
data(sir.e01, package = "fastbeta")
str(sir.e01)
plot(sir.e01)
}
