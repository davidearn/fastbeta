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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
beta <- function (t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function (t) 1e+03
mu <- function (t) 1e-03

S0 <- 5e+04
I0 <- 1e+03
constants <- c(gamma = 0.5, S0 = S0, I0 = I0, R0 = 1e+06 - S0 - I0)

n <- 2500L
X <- sir(n, beta, nu, mu, constants, stochastic = FALSE)
plot(X)

## Discard the transient
X <- window(X, start = 2250)
plot(X)

series. <- ts(cbind(unclass(X[, c("Z", "B")]), mu = mu(0)), start = 0)
series.[1L, 1L:2L] <- series.[2L, 1L:2L]
constants. <- constants[1:3]
constants.[2:3] <- floor(X[1L, c("S", "I")])
R <- fastbeta.bootstrap(10L, series., constants.)
plot(R)
plot(R, level = 0.95)
}
