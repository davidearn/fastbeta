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
\item{a}{the 0-index of the first peak in the incidence time series.}
\item{b}{the 0-index of the last peak in the incidence time series
  that is in phase with the first.}
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
\item{X}{if \code{complete = TRUE}, then a matrix with \code{b-a+1} rows
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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(graphics, pos = "package:base", verbose = FALSE)
}
beta <- function(t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function(t) 1e+03
mu <- function(t) 1e-03

S0 <- 5e+04
I0 <- 1e+03
constants <- c(gamma = 0.5, S0 = S0, I0 = I0, R0 = 1e+06 - S0 - I0)

set.seed(0)
n <- 250L
X <- sir(n, beta, nu, mu, constants)
plot(X)

series. <- ts(cbind(X[, c("Z", "B")], mu(0)), start = 0)
str(L0 <- ptpi(series., start = 4e+04, complete = FALSE))
str(L1 <- ptpi(series., start = 4e+04, complete =  TRUE))
matplot(L1$X, type = "l")
}
