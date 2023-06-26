\name{sir}
\alias{sir}
\title{Simulate Case Reports Under Forced SIR Model}
\description{
Simulates time series of the susceptible, infected, and removed
population sizes and corresponding time series of incidence and
\emph{reported} incidence.  Simulations are based on an SIR model
with user-defined forcing and a simple model for observation error.
}
\details{
\code{sir(stochastic = FALSE)} works by numerically integrating
the system of ordinary differential equations
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
between times \eqn{t = 0, \ldots, n}, where \eqn{t} is understood
to be a unitless measure of time relative to the duration of an
observation interval.  Incidence is defined as
\deqn{Z(t) = \int_{t - 1}^{t} \beta(s) S(s) I(s)\,\text{d}s}{
  Z(t) = integrate(beta(s) S(s) I(s); s = [t - 1, t])}
To compute it, the system is augmented with a fourth equation
describing \emph{cumulative} incidence
%% MJ: Why do I need \displaystyle inside of \deqn{} ... ??
\deqn{\displaystyle\frac{\text{d} c}{\text{d} t} = \beta(t) S I}{
  dc/dt = beta(t) S I}
and the \emph{augmented} system with four equations is integrated.
\code{sir(stochastic = TRUE)} works by simulating a Markov process
corresponding to the augmented system, as described in the reference.

Case reports are simulated by partitioning binomial samples
according to multinomial delays:
\deqn{Z_\text{obs}(t) = \sum\limits_{i = 0}^{d} M_{i}(t - i)}{
  Z.obs(t) = sum(M(t - i)[i]; i = 0, \ldots, d)}
where
\deqn{M(t) \sim \mathrm{Multinomial}(B(t), (q_{0}, \ldots, q_{d}))}{
  M(t) ~ Multinomial(B(t), (q[0], \ldots, q[d]))}
and
\deqn{B(t) \sim \mathrm{Binomial}(Z(t), p(t))}{
  B(t) ~ Binomial(Z(t), p(t))}
}
\usage{
sir(n, par, beta, nu, mu, stochastic = TRUE, prob = 1, delay = 1,
    useCompiled = TRUE, \dots)
}
\arguments{
\item{n}{a positive integer.  The simulation uses \code{0:n} as
  time points, yielding \code{n} observation intervals of equal
  duration.}
\item{par}{a named numeric vector of the form
  \code{c(S0 = ., I0 = ., R0 = ., gamma = .)} (or a permutation),
  specifying an initial state and a removal rate.}
\item{beta, nu, mu}{functions of one or more argument returning
  transmission, birth, and natural death rates, respectively,
  at the time points indicated by the first argument.  Arguments
  after the first must be strictly optional.}
\item{stochastic}{a logical indicating if the simulation should
  be stochastic; see \sQuote{Details}.}
\item{prob}{a numeric vector of length \code{n} such that
  \code{prob[i]} is the probability that an infection during
  observation interval \code{i} of \code{n} is eventually reported.
  \code{prob} of length 1 are recycled.}
\item{delay}{a numeric vector of positive length such that
  \code{delay[i]} is the probability that an infection is reported
  after a delay of \code{i-1} time units, given that it is eventually
  reported.  \code{delay} is passed \code{\link{sample}}, hence it
  need not sum to 1 but it must not sum to 0.}
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
\code{\link[=ts]{mts}}.  Beneath the class, it is
an \code{(n+1)}-by-\code{(4+r)} numeric matrix \code{X},
where \code{r} is 0 if both \code{prob} and \code{delay}
are missing (indicating no observation error) and 1 otherwise.

Rows correspond to time points \code{0:n}.  Columns are named
\code{S}, \code{I}, \code{R}, \code{Z}, and (if \code{r} is 1)
\code{Z.obs}.  \code{X[, 1:3]} give the state at each time, and
\code{X[, 4:5]} give the number of infections and infections
\emph{reported} during the observation interval ending at each
time.

\code{X[1, 4:5]} is \code{NA}, and \code{X[2:length(delay), 5]}
can contain incomplete information if \code{length(delay) >= 2}.
}
%\seealso{}
\references{
Cao, Y., Gillespie, D. T., & Petzold, L. R. (2007).
Adaptive explicit-implicit tau-leaping method with
automatic tau selection.
\emph{Journal of Chemical Physics},
\emph{126}(22), Article 224101, 1-9.
\doi{10.1063/1.2745299}
}
\examples{
\dontshow{ % R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
beta <- function(t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function(t) 1e+03
mu <- function(t) 1e-03

S0 <- 5e+04
I0 <- 1e+03
par <- c(S0 = S0, I0 = I0, R0 = 1e+06 - S0 - I0, gamma = 0.5)

set.seed(0)
n <- 250L
X <- sir(n, par, beta, nu, mu, stochastic = FALSE)
plot(X)

r <- 10L
Y <- replicate(r,
	sir(n, par, beta, nu, mu, prob = 0.5,
	    delay = diff(pgamma(0:8, shape = 1, scale = 1)))[, "Z.obs"])
Y <- ts(Y, start = 0)
plot(Y, plot.type = "single", col = seq_len(r))
}
