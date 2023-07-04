\name{fastbeta}
\alias{fastbeta}
\title{Estimate a time-varying infectious disease transmission rate}
\description{
Generates a discrete time estimate of the transmission rate
\eqn{\beta(t)} from an equally spaced incidence time series
and other data.
}
\details{
\code{fastbeta} works by discretizing the system of ordinary
differential equations
\deqn{
  \begin{aligned}
  \frac{\text{d} S}{\text{d} t} &= \nu(t) - \beta(t) S I - \mu(t) S \\
  \frac{\text{d} I}{\text{d} t} &= \beta(t) S I - \gamma I - \mu(t) I \\
  \frac{\text{d} R}{\text{d} t} &= \gamma I - \mu(t) R
  \end{aligned}
}{
  dS/dt = nu(t) - beta(t) S I - mu(t) S
  dI/dt = beta(t) S I - gamma I - mu(t) I
  dR/dt = gamma I - mu(t) R
}
where \eqn{t} is understood to be a unitless measure of time
relative to the duration of an observation interval, then computing
the recursion
\deqn{
  \begin{aligned}
  \beta_{t} &= \frac{Z_{t} + Z_{t + 1}}{2 S_{t} I_{t}}
  S_{t + 1} &= \frac{(1 - \frac{1}{2} \mu_{t}) S_{t} + B_{t + 1} - Z_{t + 1}}{
    1 + \frac{1}{2} \mu_{t + 1}} \\
  I_{t + 1} &= \frac{(1 - \frac{1}{2} (\gamma + \mu_{t})) I_{t} + Z_{t + 1}}{
    1 + \frac{1}{2} (\gamma + \mu_{t + 1})} \\
  \end{aligned}
}{
  beta[t] = (Z[t] + Z[t + 1]) / (2 * S[t] * I[t])
  S[t + 1] = ((1 - 0.5 * mu[t]) * S[t] + B[t + 1] - Z[t + 1]) / (1 + 0.5 * mu[t + 1])
  I[t + 1] = ((1 - 0.5 * (gamma + mu[t])) * I[t] + Z[t + 1]) / (1 + 0.5 * (gamma + mu[t + 1]))
}
where
\deqn{
  \begin{aligned}
  X_{t} &\sim X(t) \quad [X = \beta, S, I, \mu] \\
  Z_{t} &\sim \int_{t - 1}^{t} \beta(s) S(s) I(s) \, \text{d} s \\
  B_{t} &\sim \int_{t - 1}^{t} \nu(s) \, \text{d} s
  \end{aligned}
}{
  X[t] ~ X(t)    [X = beta, S, I, mu]
  Z[t] ~ integrate(beta(s) S(s) I(s), t - 1, t)
  B[t] ~ integrate(            nu(s), t - 1, t)
}
and it is understood that indexing starts at \eqn{t = 0}.
}
\usage{
fastbeta(series, constants)
}
\arguments{
\item{series}{a \dQuote{multiple time series} object, inheriting from
  class \code{\link[=ts]{mts}}, with three columns storing
  (\dQuote{parallel}, equally spaced) time series of incidence, births,
  and the per capita natural mortality rate, in that order.}
\item{constants}{a numeric vector of the form \code{c(gamma, S0, I0)},
  specifying a removal rate and an initial state, in that order.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, with three columns (named \code{beta},
\code{S}, and \code{I}) storing the result of the recursion
described in \sQuote{Details}.  It is completely parallel to
argument \code{series}, having the same \code{\link{tsp}} attribute.
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
library(stats, pos = "package:base", verbose = FALSE)
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
X <- sir(n, beta, nu, mu, constants,
         prob = 0.1, delay = diff(pgamma(0:8, 2.5)))
plot(X)

series. <- ts(cbind(X[, c("Z", "B")], mu(0)), start = 0)
constants. <- constants[1:3]
Y <- fastbeta(series., constants.)
plot(Y)
}
