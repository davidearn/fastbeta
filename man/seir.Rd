\name{seir}
\alias{seir}
\alias{sir}
\title{Simulate Infectious Disease Time Series}
\description{
Simulates incidence time series based on an \eqn{SEIR} model with
user-defined forcing and a simple model for observation error.
}
\usage{
seir(length.out = 1L,
     beta, nu, mu, sigma = gamma, gamma = 1, delta = 0,
     init, m = length(init) - n - 2L, n = 1L,
     stochastic = TRUE, prob = 1, delay = 1,
     useCompiled = TRUE, \dots)

## A basic wrapper for the 'm = 0L' case:

 sir(length.out = 1L,
     beta, nu, mu, gamma = 1, delta = 0,
     init, n = 1L,
     stochastic = TRUE, prob = 1, delay = 1,
     useCompiled = TRUE, \dots)
}
\arguments{
\item{length.out}{
  a non-negative integer indicating the time series length.}
\item{beta, nu, mu}{
  functions of one or more arguments returning transmission, birth, and
  natural death rates at the time point indicated by the first argument.
  Arguments after the first must be strictly optional.  The functions
  need not be vectorized.}
\item{sigma, gamma, delta}{
  non-negative numbers.  \code{m*sigma}, \code{n*gamma}, and
  \code{delta} are the rates of removal from each latent, infectious,
  and recovered compartment.}
\item{init}{
  a numeric vector of length \code{1+m+n+1} giving an initial state with
  compartments ordered as \eqn{(S, E, I, R)}{c(S, E, I, R)}.}
\item{m}{
  a non-negative integer indicating a number of latent stages.}
\item{n}{
  a positive integer indicating a number of infectious stages.}
\item{stochastic}{
  a logical indicating if the simulation should be stochastic; see
  \sQuote{Details}.}
\item{prob}{
  a numeric vector of length \code{n} such that \code{prob[i]} is the
  probability that an infection during interval \code{i} is eventually
  observed.  \code{prob} of length 1 is recycled.}
\item{delay}{
  a numeric vector of positive length such that \code{delay[i]} is the
  probability that an infection during interval \code{j} is observed
  during interval \code{j+i-1}, given that it is eventually observed.
  \code{delay} need not sum to 1 but must not sum to 0.}
\item{useCompiled}{
  a logical indicating if derivatives should be computed by compiled
  \proglang{C} functions rather than by \R{} functions (which \emph{may}
  be \emph{byte}-compiled).  Set to \code{FALSE} only if \code{TRUE}
  seems to cause problems, and in that case please report the problems
  with \code{\link{bug.report}(package = "fastbeta")}.}
\item{\dots}{
  optional arguments passed to \code{\link[deSolve]{lsoda}} (directly)
  or \code{\link[adaptivetau]{ssa.adaptivetau}} (via its list argument
  \code{tl.params}), depending on \code{stochastic}.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}.  Beneath the class, it is a
\code{length.out}-by-\code{1+m+n+1+2} numeric matrix with columns
\code{S}, \code{E}, \code{I}, \code{R}, \code{Z}, and \code{B}, where
\code{Z} and \code{B} specify incidence and births as the number of
infections and births since the previous time point.

If \code{prob} or \code{delay} is not missing, then there is an
additional column \code{Z.obs} specifying \emph{observed} incidence as
the number of infections observed since the previous time point.
The first \code{length(delay)} elements of this column contain partial
counts.
}
\details{
Simulations are based on an \eqn{SEIR} model with
\itemize{
\item \eqn{m} latent stages
  (\eqn{E^{i}}{E[i]}, \eqn{i = 1,\ldots,m}{for i in 1:m});
\item \eqn{n} infectious stages
  (\eqn{I^{j}}{I[j]}, \eqn{j = 1,\ldots,n}{for j in 1:n});
\item time-varying rates \eqn{\beta}, \eqn{\nu}, and \eqn{\mu} of
  transmission, birth, and natural death; and
\item constant rates \eqn{m \sigma}, \eqn{n \gamma}, and \eqn{\delta} of
  removal from each latent, infectious, and recovered compartment, where
  removal from the recovered compartment implies return to the
  susceptible compartment (loss of immunity).
}
\code{seir(stochastic = FALSE)} works by numerically integrating the
system of ordinary differential equations
\deqn{
  \begin{alignedat}{10}
  \text{d} & S         &{} / \text{d} t
  &{} = {}&     \delta &R     &{} - ( && \lambda(t) &{} + \mu(t)) S         &{} + \nu(t) \\
  \text{d} & E^{    1} &{} / \text{d} t
  &{} = {}& \lambda(t) &S     &{} - ( &&   m \sigma &{} + \mu(t)) E^{    1} &{}          \\
  \text{d} & E^{i + 1} &{} / \text{d} t
  &{} = {}&   m \sigma &E^{i} &{} - ( &&   m \sigma &{} + \mu(t)) E^{i + 1} &{}          \\
  \text{d} & I^{    1} &{} / \text{d} t
  &{} = {}&   m \sigma &E^{m} &{} - ( &&   n \gamma &{} + \mu(t)) I^{    1} &{}          \\
  \text{d} & I^{j + 1} &{} / \text{d} t
  &{} = {}&   n \gamma &I^{j} &{} - ( &&   n \gamma &{} + \mu(t)) I^{j + 1} &{}          \\
  \text{d} & R         &{} / \text{d} t
  &{} = {}&   n \gamma &I^{n} &{} - ( &&     \delta &{} + \mu(t)) R         &{}
  \end{alignedat} \\
  \lambda(t) = \beta(t) \sum_{j} I^{j}
}{
  dS       /dt =     delta * R    - ( lambda(t) + mu(t) ) * S        + nu(t)
  dE[    1]/dt = lambda(t) * S    - ( m * sigma + mu(t) ) * E[    1]
  dE[i + 1]/dt = m * sigma * E[i] - ( m * sigma + mu(t) ) * E[i + 1]
  dI[    1]/dt = m * sigma * E[m] - ( n * gamma + mu(t) ) * I[    1]
  dI[j + 1]/dt = n * gamma * I[j] - ( n * gamma + mu(t) ) * I[j + 1]
  dR       /dt = n * gamma * I[n] - (     delta + mu(t) ) * R

  lambda(t) = beta(t) * sum(I)
}
where it is understood that the independent variable \eqn{t} is a
unitless measure of time relative to an observation interval.  To get
time series of incidence and births, the system is augmented with two
equations describing \emph{cumulative} incidence and births
\deqn{
  \begin{aligned}
  \text{d} Z / \text{dt} &{} = \lambda(t) S
  \text{d} B / \text{dt} &{} = \nu(t)
  \end{aligned}
}{
  dZ/dt = lambda(t) * S
  dB/dt = nu(t)
}
and the \emph{augmented} system is numerically integrated.
Observed incidence is simulated from incidence by scaling the latter
by \code{prob} and convolving the result with \code{delay}.

\code{seir(stochastic = TRUE)} works by simulating a Markov process
corresponding to the augmented system, as described in the reference.
Observed incidence is simulated from incidence by binning binomial
samples taken with probabilities \code{prob} over future observation
intervals according to multinomial samples taken with probabilities
\code{delay}.
}
\seealso{\code{\link{seir.library}}.}
\references{
Cao, Y., Gillespie, D. T., & Petzold, L. R. (2007).
Adaptive explicit-implicit tau-leaping method with
automatic tau selection.
\emph{Journal of Chemical Physics},
\emph{126}(22), Article 224101, 1-9.
\doi{10.1063/1.2745299}
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
beta <- function (t, a = 1e-01, b = 1e-05) b * (1 + a * sinpi(t / 26))
nu   <- function (t) 1e+03
mu   <- function (t) 1e-03

sigma <- 0.5
gamma <- 0.5
delta <- 0

init <- c(S = 50200, E = 1895, I = 1892, R = 946011)

length.out <- 250L
prob <- 0.1
delay <- diff(pgamma(0:8, 2.5))

set.seed(0L)
X <- seir(length.out, beta, nu, mu, sigma, gamma, delta, init,
          prob = prob, delay = delay, epsilon = 0.002)
##                                    ^^^^^^^
## default epsilon = 0.05 allows too big leaps => spurious noise
##
str(X)
plot(X)

r <- 10L
Y <- do.call(cbind, replicate(r, simplify = FALSE,
	seir(length.out, beta, nu, mu, sigma, gamma, delta, init,
	     prob = prob, delay = delay, epsilon = 0.002)[, "Z.obs"]))
str(Y) # FIXME: stats:::cbind.ts mangles dimnames
plot(window(Y, start = tsp(Y)[1L] + length(delay) / tsp(Y)[3L]),
     ##        ^^^^^
     ## discards points showing edge effects due to 'delay'
     ##
     plot.type = "single", col = seq_len(r), ylab = "case reports")
}
