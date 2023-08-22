\name{fastbeta.bootstrap}
\alias{fastbeta.bootstrap}
\title{Parametric bootstrapping}
\description{
A simple wrapper around \code{\link{fastbeta}}, using it to generate
a \dQuote{primary} estimate of a transmission rate \eqn{\beta(t)}
and \code{r} bootstrap estimates.  Bootstrap estimates are computed
for incidence time series simulated using \code{\link{sir}}, with
transmission rate defined as the linear interpolant of the primary
estimate.
}
%\details{}
\usage{
fastbeta.bootstrap(r, series, constants, \dots)
}
\arguments{
\item{r}{the number of replications.}
\item{series, constants}{see \code{\link{fastbeta}}.}
\item{\dots}{optional arguments passed to \code{\link{sir}}
  and/or \code{\link{deconvolve}}.  Both take optional arguments
  \code{prob} and \code{delay}.  When \code{prob} is supplied
  but not \code{delay}, \code{sir} and \code{deconvolve} receive
  \code{prob} as is.  When both are supplied, \code{\link{sir}}
  receives \code{prob} as is, whereas \code{\link{deconvolve}}
  receives \code{prob} augmented with \code{length(delay)-1}
  ones.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, with \code{r+1} columns storing the primary
and bootstrap estimates, in that order.
It is completely parallel to argument \code{series}, having the same
\code{\link{tsp}} attribute.
}
%\seealso{}
%\references{}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(graphics, pos = "package:base", verbose = FALSE)
library(   utils, pos = "package:base", verbose = FALSE)
}
data(sir.e01, package = "fastbeta")
a <- attributes(sir.e01)
str(sir.e01)
plot(sir.e01)

## Now suppose that you have perfect knowledge except for incidence,
## which you observe imperfectly
series <- cbind(sir.e01[, c("Z.obs", "B")], mu = a[["mu"]](0))
constants <- c(gamma = a[["gamma"]],
               S0 = sir.e01[[1L, "S"]],
               I0 = sir.e01[[1L, "I"]])

R <- fastbeta.bootstrap(10L, series, constants)
str(R)
plot(R)
plot(R, level = 0.95)
}
