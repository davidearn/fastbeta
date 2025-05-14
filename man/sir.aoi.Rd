\name{sir.aoi}
\alias{sir.aoi}
\alias{summary.sir.aoi}
\title{Solve the SIR Equations Structured by Age of Infection}
\description{
Numerically integrates the SIR equations with rates of transmission and
recovery structured by age of infection.
}
\usage{
sir.aoi(from = 0, to = from + 1, by = 1,
        R0, ell = 1, n = max(length(R0), length(ell)),
        init = c(1 - init.infected, init.infected),
        init.infected = .Machine[["double.neg.eps"]],
        weights = rep(c(1, 0), c(1L, n - 1L)),
        H = identity, Hargs = list(), root = NULL,
        aggregate = FALSE, \dots)

\method{summary}{sir.aoi}(object, name = "Y", tol = 1e-6, \dots)
}
\arguments{
\item{from, to, by}{
  passed to \code{\link{seq.int}} in order to generate an increasing,
  equally spaced vector of time points in units of the mean time spent
  infectious.}
\item{R0}{
  a numeric vector of length \code{n} such that \code{sum(R0)} is the
  basic reproduction number and \code{R0[j]} is the contribution of
  infected compartment \code{j}.  Otherwise, a numeric vector of length
  1, handled as equivalent to \code{rep(R0/n, n)}.}
\item{ell}{
  a numeric vector of length \code{n} such that \code{ell[j]} is the
  ratio of the mean time spent in infected compartment \code{j} and the
  mean time spent infectious; internally, \code{ell/sum(ell[R0 > 0])}
  is used, hence \code{ell} is determined only up to a positive factor.
  Otherwise (and by default), a numeric vector of length 1, handled as
  equivalent to \code{rep(1, n)}.}
\item{n}{
  a positive integer giving the number of infected compartments.
  Setting \code{n} and thus overriding the default expression is
  necessary only if the lengths of \code{R0} and \code{ell} are both 1.}
\item{init}{
  a numeric vector of length 2 giving initial susceptible and infected
  proportions.}
\item{init.infected}{
  a number in \eqn{(0, 1]} used only to define the default expression
  for \code{init}; see \sQuote{Usage}.}
\item{weights}{
  a numeric vector of length \code{n} containing non-negative weights,
  defining the initial distribution of infected individuals among the
  infected compartments.  By default, all infected individuals occupy
  the first compartment.}
\item{H}{
  a function returning a numeric vector of length equal to that of the
  first formal argument.  The body must be symbolically differentiable
  with respect to the first formal argument; see \code{\link{D}}.}
\item{Hargs}{
  a list of arguments passed to \code{H}.  In typical usage,
  \code{H} defines a parametric family of functions of one variable,
  and \code{Hargs} supplies parameter values.  For example:
  \code{H = function(x, h) x^h}, \code{Hargs = list(h = 0.996)}.}
\item{root}{
  a function returning a numeric vector of length 1, with formal
  arguments \code{(tau, S, I, Y, dS, dI, dY, R0, ell)} (or a subset);
  otherwise, \code{NULL}.}
\item{aggregate}{
  a logical indicating if infected compartments should be aggregated.}
\item{\dots}{
  optional arguments passed to \code{\link[deSolve]{lsoda}}; use
  \code{atol}, \code{rtol}, and \code{hmax} for control of local error.}
\item{object}{
  an \R{} object inheriting from class \code{sir.aoi}, typically
  the value of a call to \code{sir.aoi}.}
\item{name}{
  a character string, either \code{"Y"} or \code{"I"}.  Tail exponents
  of \eqn{V} are approximated, where \eqn{V = Y} if \code{name = "Y"}
  and \eqn{V = \sum_{j} I_{j}}{V = sum(I)} (prevalence) if
  \code{name = "I"}.}
\item{tol}{
  a positive number giving an upper bound on the relative change (from
  one time point to the next) in the slope of \eqn{\log(V)}, defining
  time windows in which growth or decay of \eqn{V} is considered to be
  exponential.}
}
\value{
\item{root = NULL}{
  a \dQuote{multiple time series} object, inheriting from class
  \code{sir.aoi} and transitively from class \code{\link[=ts]{mts}}.
  Beneath the class, it is a
  \code{length(seq(from, to, by))}-by-\code{(1+n+1)} numeric matrix
  of the form \code{cbind(S, I, Y)}.}
\item{root = function (tau, S, I, Y, dS, dI, dY, R0, ell)}{
  a numeric vector of length \code{1+1+n+1} of the form
  \code{c(tau, S, I, Y)} storing the root of the function \code{root} in
  units of the mean time spent infectious and the state at that time.
  Attribute \code{curvature} stores the curvature of \eqn{Y} at the
  root.  If a root is not found between times \code{from} and \code{to},
  then the value is \code{NULL}.}

If \code{aggregate = TRUE}, then infected compartments are aggregated so
that the number of columns (elements, if \code{root} is a function)
named \code{I} is 1 rather than \code{n}.  This column or element stores
prevalence, the proportion of the population that is infected.
For convenience, there are 4 additional columns (elements) named
\code{I.E}, \code{I.I}, \code{foi}, and \code{inc}.  These store the
non-infectious and infectious components of prevalence (so that
\code{I.E + I.I = I}), the force of infection, and incidence.

The method for \code{\link{summary}} returns a numeric vector of length
2 containing the left and right \dQuote{tail exponents} of the variable
\eqn{V} indicated by \code{name}, defined as the asymptotic values of
the slope of \eqn{\log(V)}.  \code{NaN} elements indicate that a tail
exponent cannot be approximated because the time window represented by
\code{object} does not cover enough of the tail, where the meaning of
\dQuote{enough} is set by \code{tol}.
}
\details{
The SIR equations with rates of transmission and recovery structured by
age of infection are
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} t
  &{} = -(\textstyle\sum_{j} \beta_{j} I_{j}) H(S) \\
  \text{d} & I_{    1} &{} / \text{d} t
  &{} = (\textstyle\sum_{j} \beta_{j} I_{j}) H(S) - \gamma_{1} I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} t
  &{} = \gamma_{j} I_{j} - \gamma_{j + 1} I_{j + 1} \\
  \text{d} & R         &{} / \text{d} t
  &{} = \gamma_{n} I_{n}
  \end{alignedat}
}{
  dS       /dt = -sum(beta * I) * H(S)
  dI[    1]/dt = sum(beta * I) * H(S) - gamma[1] * I[1]
  dI[j + 1]/dt = gamma[j] * I[j] - gamma[j + 1] * I[j + 1]
  dR       /dt = gamma[n] * I[n]
}
where \eqn{S, I_{j}, R \ge 0}{min(S, I, R) >= 0},
\eqn{S + \sum_{j} I_{j} + R = 1}{sum(S, I, R) = 1}, and \eqn{H} is the
susceptible heterogeneity function.  In general, \eqn{H} is nonlinear;
in the standard SIR equations, \eqn{H} is the identity function.

Nondimensionalization using parameters
\eqn{\mathcal{R}_{0,j} = \beta_{j}/\gamma_{j}}{R0 = beta/gamma} and
\eqn{\ell_{j} = (1/\gamma_{j})/\sum_{j:\mathcal{R}_{0,j} > 0} (1/\gamma_{j})}{ell = (1/gamma)/sum(1/gamma[R0 > 0])}
and time unit
\eqn{\tau = t/\sum_{j:\mathcal{R}_{0,j} > 0} (1/\gamma_{j})}{tau = t/sum(1/gamma[R0 > 0])},
gives
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} \tau
  &{} = -(\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) I_{j}) H(S) \\
  \text{d} & I_{    1} &{} / \text{d} \tau
  &{} = (\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) I_{j}) H(S) - (1/\ell_{1}) I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} \tau
  &{} = (1/\ell_{j}) I_{j} - (1/\ell_{j+1}) I_{j + 1} \\
  \text{d} & R         &{} / \text{d} \tau
  &{} = (1/\ell_{n}) I_{n} \\
  \end{alignedat}
}{
  dS       /dtau = -sum((R0/ell) * I) * H(S)
  dI[    1]/dtau = sum((R0/ell) * I) * H(S) - (1/ell[1]) * I[1]
  dI[j + 1]/dtau = (1/ell[j]) * I[j] - (1/ell[j + 1]) * I[j + 1]
  dR       /dtau = (1/ell[n]) * I[n]
}
\code{sir.aoi} works with the nondimensional equations, dropping the
last equation (which is redundant given
\eqn{R = 1 - S - \sum_{j} I_{j}}{R = 1 - sum(S, I)}) and augments the
resulting system of \eqn{1 + n} equations with a new equation
\deqn{
  \text{d}Y/\text{d}\tau = (\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) I_{j}) (H(S) - H(1/\textstyle\sum_{j} \mathcal{R}_{0, j}))
}{
  dY/tau = sum((R0/ell) * I) * (H(S) - H(1/sum(R0)))
}
due to the usefulness of the solution \eqn{Y} in applications.
}
\note{
\code{sir.aoi} is not a special case of \code{\link{sir}} nor a
generalization.  The two functions were developed independently and for
different purposes: \code{sir.aoi} to validate analytical results
concerning the SIR equations as formulated here, \code{sir} to simulate
incidence time series suitable for testing \code{\link{fastbeta}}.
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(   stats, pos = "package:base", verbose = FALSE)
library(graphics, pos = "package:base", verbose = FALSE)
}
if (requireNamespace("deSolve")) withAutoprint({

to <- 100; by <- 0.01; R0 <- c(0, 16); ell <- c(0.5, 1)

peak <- sir.aoi(to = to, by = by, R0 = R0, ell = ell,
                root = function (S, R0) sum(R0) * S - 1,
                aggregate = TRUE)
peak

to <- 4 * peak[["tau"]] # a more principled endpoint

soln <- sir.aoi(to = to, by = by, R0 = R0, ell = ell,
                aggregate = TRUE, atol = 0, rtol = 1e-12)
##                                ^^^^      ^^^^
## 'atol', 'rtol', 'hmax' are passed to deSolve::lsoda

head(soln)
tail(soln)

plot(soln) # dispatching stats:::plot.ts
plot(soln, log = "y")

(lamb <- summary(soln)) # left and right tail exponents

xoff <- function (x, k) x - x[k]
k <- c(16L, nrow(soln)) # c(1L, nrow(soln)) worse due to transient

plot(soln[, "Y"], log = "y", ylab = "Y")
abline(v = peak[["tau"]], h = peak[["Y"]],
       lty = 2, lwd = 2, col = "red")
for (i in 1:2)
	lines(soln[k[i], "Y"] * exp(lamb[i] * xoff(time(soln), k[i])),
	      lty = 2, lwd = 2, col = "red")

wrap <-
function (root)
	sir.aoi(to = to, by = by, R0 = R0, ell = ell,
	        root = root, aggregate = TRUE)
Ymax <- peak[["Y"]]

## NB: want *simple* roots, not *multiple* roots
F <- list(function (Y) (Y - Ymax * 0.5)  ,
          function (Y) (Y - Ymax * 0.5)^2,
          function (Y) (Y - Ymax      )  ,
          function (Y) (Y - Ymax      )^2)
lapply(F, wrap)

## NB: critical values can be attained twice
F <- list(function (Y, dY) if (dY > 0) Y - Ymax * 0.5 else 1,
          function (Y, dY) if (dY < 0) Y - Ymax * 0.5 else 1)
lapply(F, wrap)

})
}
