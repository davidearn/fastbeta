\name{ptpi}
\alias{ptpi}
\title{Peak to peak iteration}
\description{
Estimates the initial susceptible population size from a periodic,
equally spaced incidence time series and other data.
}
\details{
\code{ptpi} works by computing the recursion
\deqn{
  \begin{aligned}
  S_{t + 1}^{(i)} &= \frac{(1 - \frac{1}{2} \mu_{t}) S_{t}^{(i)} + B_{t + 1} - Z_{t + 1}}{1 + \frac{1}{2} \mu_{t + 1}} \quad [t = a, \ldots, b - 1] \\
  S_{a}^{(i + 1)} &= S_{b}^{(i)}
  \end{aligned}
}{
  S[t + 1, i] = ((1 - 0.5 * mu[t]) * S[t, i] + B[t + 1] - Z[t + 1]) / (1 + 0.5 * mu[t + 1])    [t = a, \ldots, b - 1]
  S[a, i + 1] = S[b, i]
}
starting with \eqn{S_{a}^{(0)}}{S[a, i]} equal to \code{start} and
ending when
\deqn{
  \left|\frac{S_{b}^{(i)} - S_{a}^{(i)}}{S_{a}^{(i)}}\right| < \mathtt{tol}
}{
  abs((S[b, i] - S[a, i]) / S[a, i]) < tol
}
or when \eqn{i} reaches \code{iter.max}.
The estimated value of the initial susceptible population size is
\eqn{S_{0}^{(i + 1)}}{S[0, i + 1]}, which is back-calculated from
\eqn{S_{a}{(i + 1)}}{S[a, i + 1]} by reversing the recursion.
}
\usage{
ptpi(series, a = 0L, b = nrow(series) - 1L,
     start, tol = 1e-03, iter.max = 32L, complete = FALSE, \dots)
}
\arguments{
\item{series}{a \dQuote{multiple time series} object, inheriting from
  class \code{\link[=ts]{mts}}, with three columns storing
  (\dQuote{parallel}, equally spaced) time series of incidence, births,
  and the per capita natural mortality rate, in that order.}
\item{a}{the time of the first peak in the incidence time series.
  It is rounded internally to generate a 0-index of rows of \code{series}.}
\item{b}{the time of the last peak in the incidence time series
  that is in phase with the first.
  It is rounded internally to generate a 0-index of rows of \code{series}.}
\item{start}{a starting value for the number of susceptible individuals
  at index \code{a}.}
\item{tol}{a tolerance indicating a stopping condition;
  see \sQuote{Details}.}
\item{iter.max}{the maximum number of iterations.}
\item{complete}{a logical flag indicating if the result should
  preserve the number of susceptible individuals from index \code{a}
  to index \code{b} in each iteration.}
\item{\dots}{optional arguments passed to \code{\link{deconvolve}},
  if the first column of \code{series} represents \emph{reported}
  incidence or mortality rather than incidence.}
}
\value{
A list with elements:
\item{value}{the estimated value of the initial susceptible population size.}
\item{delta}{the relative difference computed in the last iteration.}
\item{iter}{the number of iterations performed.}
\item{X}{if \code{complete = TRUE}, then a \dQuote{multiple time series}
  object, inheriting from class \code{\link[=ts]{mts}}, with \code{b-a+1} rows
  and \code{iter} columns storing the value of \eqn{S_{t}^{(i)}}{S[t, i]}
  for \eqn{t = a, \ldots, b}
  and \eqn{i = 0, \ldots, \mathtt{iter} - 1}{i = 0, \ldots, iter - 1}.}
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

series <- cbind(sir.e01[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu") # FIXME: stats:::cbind.ts mangles dimnames

plot(series[, "Z"])
a <- 9; b <- 218
abline(v = c(a, b), lty = 2)

L <- ptpi(series, a = a, b = b, start = 4e+04, complete = TRUE)
str(L)
plot(L$X, plot.type = "single")
lines(sir.e01[, "S"], col = "red", lwd = 4)
abline(h = L$value, col = "blue", lwd = 4)

## The relative error:
L$value / sir.e01[[1L, "S"]] - 1
}
