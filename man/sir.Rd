\name{sir}
\alias{sir}
\title{Simulate infectious disease time series}
\description{
Simulates time series of the susceptible, infected, and removed
population sizes and corresponding time series of births, incidence,
and \emph{reported} incidence.  Simulations are based on an SIR model
with user-defined forcing and a simple model for observation error.
}
\details{
\code{sir(stochastic = FALSE)} works by numerically integrating the
system of ordinary differential equations
\deqn{
  \begin{aligned}
  \frac{\text{d} S}{\text{d} t} &= \nu(t) - \beta(t) S I + \delta R - \mu(t) S \\
  \frac{\text{d} I}{\text{d} t} &= \beta(t) S I - \gamma I - \mu(t) I \\
  \frac{\text{d} R}{\text{d} t} &= \gamma I - \delta R - \mu(t) R
  \end{aligned}
}{
  dS/dt = nu(t) - beta(t) S I + delta R - mu(t) S
  dI/dt = beta(t) S I - gamma I - mu(t) I
  dR/dt = gamma I - delta R - mu(t) R
}
between times \eqn{t = 0, \ldots, n}, where \eqn{t} is understood
to be a unitless measure of time relative to the duration of an
observation interval.  To generate time series of births and incidence
\deqn{
  \begin{aligned}
  B(t) &= \int_{t - 1}^{t} \nu(s)\,\text{d}s \\
  Z(t) &= \int_{t - 1}^{t} \beta(s) S(s) I(s)\,\text{d}s
  \end{aligned}
}{
  B(t) = integrate(            nu(s), t - 1, t)
  Z(t) = integrate(beta(s) S(s) I(s), t - 1, t)
}
the system is augmented with two additional equations describing
\emph{cumulative} incidence and \emph{cumulative} births (with right
hand sides given by the integrands above), and the \emph{augmented}
system with five equations is integrated.  Case reports are simulated
by scaling incidence by \code{prob} and convolving the result with
\code{delay}.

\code{sir(stochastic = TRUE)} works by simulating a Markov process
corresponding to the augmented system, as described in the reference.
Case reports are simulated from incidence by binning binomial samples
taken with probabilities \code{prob} over future observation intervals
according to multinomial samples taken with probabilities \code{delay}.
}
\usage{
sir(n, beta, nu, mu, constants, stochastic = TRUE,
    prob = 1, delay = 1, useCompiled = TRUE, \dots)
}
\arguments{
\item{n}{a positive integer.  The simulation uses \code{0:n} as
  time points, yielding \code{n} observation intervals of equal
  duration.}
\item{beta, nu, mu}{functions of one or more argument returning
  transmission, birth, and natural death rates, respectively,
  at the time point indicated by the first argument.  Arguments
  after the first must be strictly optional.  The functions need
  not be vectorized.}
\item{constants}{a numeric vector of the form
  \code{c(S0, I0, R0, gamma, delta)}, specifying an initial state
  and rates of removal and loss of immunity, in that order.}
\item{stochastic}{a logical indicating if the simulation should
  be stochastic; see \sQuote{Details}.}
\item{prob}{a numeric vector of length \code{n} such that
  \code{prob[i]} is the probability that an infection during
  interval \code{i} is eventually reported.  \code{prob} of
  length 1 is recycled.}
\item{delay}{a numeric vector of positive length such that
  \code{delay[i]} is the probability that an infection during
  interval \code{j} is reported during interval \code{j+i-1},
  given that it is eventually reported.  \code{delay} need
  not sum to 1 but must not sum to 0.}
\item{useCompiled}{a logical indicating if derivatives should be
  computed by compiled \proglang{C} functions rather than by \R{}
  functions (which \emph{may} be \emph{byte}-compiled).  Set to
  \code{FALSE} only if \code{TRUE} seems to cause problems, and
  in that case please report the problems with
  \code{\link{bug.report}(package="fastbeta")}.}
\item{\dots}{optional arguments passed to
  \code{\link[deSolve]{ode}} (directly) or
  \code{\link[adaptivetau]{ssa.adaptivetau}}
  (via its list argument \code{tl.params}),
  depending on \code{stochastic}.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}.  Beneath the class,
it is an \code{(n+1)}-by-\code{(5+r)} numeric matrix \code{X},
where \code{r} is 0 if both \code{prob} and \code{delay}
are missing (indicating no observation error) and 1 otherwise.

Rows correspond to time points \code{0:n}.  Columns are named
\code{S}, \code{I}, \code{R}, \code{B}, \code{Z}, and (if \code{r}
is 1) \code{Z.obs}.  \code{X[, 1:3]} give the state at each time,
and \code{X[, 4:6]} give the number of births, infections, and
infections \emph{reported} during the observation interval ending
at each time.

\code{X[1, 4:6]} is \code{NA}, and \code{X[2:length(delay), 6]}
can contain incomplete information if \code{length(delay) >= 2}.
}
\seealso{\code{\link{sir.library}}.}
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
beta <- function (t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function (t) 1e+03
mu <- function (t) 1e-03

S0 <- 5e+04
I0 <- 1e+03
R0 <- 1e+06 - S0 - I0
constants <- c(S0 = S0, I0 = I0, R0 = R0, gamma = 0.5, delta = 0)

n <- 250L
prob <- 0.1
delay <- diff(pgamma(0:8, 2.5))

set.seed(0L)
X <- sir(n, beta, nu, mu, constants, prob = prob, delay = delay)
str(X)
plot(X)

r <- 10L
Y <- do.call(cbind, replicate(r, simplify = FALSE,
	sir(n, beta, nu, mu, constants, prob = prob, delay = delay)[, "Z.obs"]))
str(Y) # FIXME: stats:::cbind.ts mangles dimnames
plot(Y, plot.type = "single", col = seq_len(r), ylab = "Case reports")
}
