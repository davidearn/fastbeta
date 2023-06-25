\name{sir}
\alias{sir}
\title{Simulate Case Reports Under Forced SIR Model}
\description{
  TODO
}
\details{
  It is understood that time is everywhere measured in units of
  the observation interval duration.
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
    observation interval \code{i} is eventually reported.}
  \item{delay}{a numeric vector of positive length \code{d} such that
    \code{delay[d-i]} is the probability that an infection is reported
    after \code{i} time units, given that it is eventually reported.
    \code{delay} is passed \code{\link{sample}}, hence it need not sum
    to 1 but must not sum to 0.}
  \item{useCompiled}{a logical indicating if derivatives should be
    computed by compiled \proglang{C} functions rather than by \R{}
    functions (which \emph{may} be \emph{byte}-compiled).  Set to
    \code{FALSE} only if \code{TRUE} seems to cause problems, and in
    that case please report the problems with
    \code{\link{bug.report}(package="fastbeta")}.}
  \item{\dots}{optional arguments passed to
    \code{\link[adaptivetau]{ssa.adaptivetau}}
    (via its list argument \code{tl.params}) or
    \code{\link[deSolve]{ode}} (directly),
    depending on \code{stochastic}.}
}
\value{
  A \dQuote{multiple time series} object, inheriting from
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

  \code{X[1, 4:5]} is \code{NA}, and \code{X[2:d, 5]} can contain
  incomplete information if \code{d = length(delay) >= 2}.
}
%\seealso{}
%\references{}
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
