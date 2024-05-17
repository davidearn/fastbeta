\name{seir.library}
\alias{seir.library}
\alias{seir.ts01}
\alias{seir.ts02}
\docType{data}
\title{Often Used Simulations}
\description{
Infectious disease time series simulated using \code{\link{seir}}, for
use primarily in examples, tests, and vignettes.  Users should not rely
on simulation details, which may change between package versions.
}
\usage{
data(seir.ts01, package = "fastbeta")
data(seir.ts02, package = "fastbeta")
}
\format{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, always a subset of the result of a call to
\code{\link{seir}}, discarding transient behaviour.  Simulation
parameters may be preserved as attributes.
}
\source{
Scripts sourced by \code{\link{data}} to reproduce the simulations are
located in subdirectory \file{data} of the \pkg{fastbeta} installation;
see, e.g.
\code{\link{system.file}("data", "seir.ts01.R", package = "fastbeta")}.
}
\seealso{\code{\link{seir}}.}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
data(seir.ts01, package = "fastbeta")
str(seir.ts01)
plot(seir.ts01)

data(seir.ts02, package = "fastbeta")
str(seir.ts02)
plot(seir.ts02)
}
