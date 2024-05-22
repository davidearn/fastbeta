\name{fastbeta.bootstrap}
\alias{fastbeta.bootstrap}
\title{Parametric Bootstrapping}
\description{
A simple wrapper around \code{\link{fastbeta}} using it to generate a
\dQuote{primary} estimate of a time-varying transmission rate and
\code{r} bootstrap estimates.  Bootstrap estimates are computed for
incidence time series simulated using \code{\link{seir}}, with
transmission rate defined as the linear interpolant of the primary
estimate.
}
\usage{
fastbeta.bootstrap(r,
                   series, sigma = gamma, gamma = 1, delta = 0,
                   m = length(init) - n - 2L, n = 1L, init, \dots)
}
\arguments{
\item{r}{
  a non-negative integer indicating a number of replications.}
\item{series}{
  a \dQuote{multiple time series} object, inheriting from class
  \code{\link[=ts]{mts}}, with three columns storing (\dQuote{parallel},
  equally spaced) time series of incidence, births, and the per capita
  natural mortality rate, in that order.}
\item{sigma, gamma, delta}{
  non-negative numbers.  \code{m*sigma}, \code{n*gamma}, and
  \code{delta} are the rates of removal from each latent, infectious,
  and recovered compartment.}
\item{m}{
  a non-negative integer indicating a number of latent stages.}
\item{n}{
  a positive integer indicating a number of infectious stages.}
\item{init}{
  a numeric vector of length \code{1+m+n+1} giving an initial state with
  compartments ordered as \eqn{(S, E, I, R)}{c(S, E, I, R)}.}
\item{\dots}{
  optional arguments passed to \code{\link{seir}} and/or
  \code{\link{deconvolve}}.  Both take optional arguments \code{prob}
  and \code{delay}.  When \code{prob} is supplied but not \code{delay},
  \code{seir} and \code{deconvolve} receive \code{prob} as is.  When
  both are supplied, \code{\link{seir}} receives \code{prob} as is,
  whereas \code{\link{deconvolve}} receives \code{prob} augmented with
  \code{length(delay)-1} ones.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, with \code{1+r} columns storing the one primary
and \code{r} bootstrap estimates.  It is completely parallel to argument
\code{series}, having the same \code{\link{tsp}} attribute.
}
%\details{}
%\seealso{}
%\references{}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
if (requireNamespace("adaptivetau")) withAutoprint({

data(seir.ts02, package = "fastbeta")
a <- attributes(seir.ts02)
str(seir.ts02)
plot(seir.ts02)

## We suppose that we have perfect knowledge of incidence,
## births, and the data-generating parameters
series <- cbind(seir.ts02[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu") # FIXME: stats:::cbind.ts mangles dimnames

args <- c(list(r = 100L, series = series),
          a[c("sigma", "gamma", "delta", "m", "n", "init")])
str(args)

R <- do.call(fastbeta.bootstrap, args)
str(R)
plot(R)
plot(R, level = 0.95)

})
}
