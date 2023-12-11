\name{ptpi}
\alias{ptpi}
\title{Peak to Peak Iteration}
\description{
Estimates the initial sizes of the susceptible, infected, and
removed populations from a periodic, equally spaced incidence
time series and other data.  Interpret with care, notably when
supplying time series that are only \dQuote{roughly} periodic.
}
\details{
\code{ptpi} works by computing the iteration described in
\code{\link{fastbeta}} from time \eqn{t = a} to time \eqn{t = b},
iteratively, until the relative difference between the states
at times \eqn{a} and \eqn{b} descends below a tolerance.
The state at time \eqn{a} in the first iteration is specified
by the user.  In subsequent iterations, it is the value of
the state at time \eqn{b} calculated in the previous iteration.

If \code{backcalc = FALSE}, then \code{ptpi} returns a list
with component \code{value} equal to the state at time \eqn{b}
in the last iteration.  By periodicity, this value estimates
the \dQuote{true} state at time \eqn{a}.

If \code{backcalc = TRUE}, then \code{value} is back-calculated
so that it estimates the \dQuote{true} state at time \eqn{0}.
This works by inverting the transformation defining one step of
the iteration, hence the back-calculated result can be nonsense
if the inverse problem is ill-conditioned.
}
\usage{
ptpi(series, constants, a = 0L, b = nrow(series) - 1L,
     tol = 1e-03, iter.max = 32L,
     complete = FALSE, backcalc = FALSE, \dots)
}
\arguments{
\item{series}{a \dQuote{multiple time series} object, inheriting from
  class \code{\link[=ts]{mts}}, with three columns storing
  (\dQuote{parallel}, equally spaced) time series of incidence, births,
  and the per capita natural mortality rate, in that order.}
\item{constants}{a numeric vector of the form
  \code{c(Sa, Ia, Ra, gamma, delta)}, specifying a starting value for
  the state at time \code{a} and rates of removal and loss of immunity,
  in that order.}
\item{a}{the time of the first peak in the incidence time series.
  It is rounded internally to generate a 0-index of rows of \code{series}.}
\item{b}{the time of the last peak in the incidence time series
  that is in phase with the first.
  It is rounded internally to generate a 0-index of rows of \code{series}.}
\item{tol}{a tolerance indicating a stopping condition;
  see \sQuote{Details}.}
\item{iter.max}{the maximum number of iterations.}
\item{complete}{a logical indicating if the result should preserve
  the state at times \code{a, \ldots, b} in each iteration.}
\item{backcalc}{a logical indicating if the state at time \code{0}
  should be back-calculated.}
\item{\dots}{optional arguments passed to \code{\link{deconvolve}},
  if the first column of \code{series} represents \emph{reported}
  incidence or mortality rather than incidence.}
}
\value{
A list with elements:
\item{value}{the estimated value of the initial state, which is the
  state at time \code{a} if \code{backcalc = FALSE} and the state at
  time \code{0} if \code{backcalc = TRUE}.}
\item{delta}{the relative difference computed in the last iteration.}
\item{iter}{the number of iterations performed.}
\item{x}{if \code{complete = TRUE},
  then a \dQuote{multiple time series} object, inheriting from class
  \code{\link[=ts]{mts}}, with dimensions \code{c(b-a+1, 3, iter)}.
  \code{x[, , i]} gives the state at times \code{a, \ldots, b} in
  iteration \code{i}.}
}
%\seealso{}
\references{
Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
Fast estimation of time-varying infectious disease transmission rates.
\emph{PLOS Computational Biology},
\emph{16}(9), Article e1008124, 1-39.
\doi{10.1371/journal.pcbi.1008124}
}
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

## We suppose that we have perfect knowledge of incidence,
## births, and the data-generating parameters, except for
## the initial state, which we "guess"
series <- cbind(sir.e01[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu") # FIXME: stats:::cbind.ts mangles dimnames
constants <- c(Sa = sir.e01[[1L, "S"]],
               Ia = sir.e01[[1L, "I"]],
               Ra = sir.e01[[1L, "R"]],
               gamma = a[["gamma"]],
               delta = a[["delta"]])

plot(series[, "Z"])
a <- 9; b <- 218
abline(v = c(a, b), lty = 2)

L <- ptpi(series, constants, a = a, b = b, complete = TRUE)
str(L)
plot(L[["S"]], plot.type = "single")
lines(sir.e01[, "S"], col = "red", lwd = 4)
abline(h = L$value[["S"]], col = "blue", lwd = 4)

## The relative error:
L$value[["S"]] / sir.e01[[1L, "S"]] - 1
}
