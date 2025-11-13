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
        R0, ell = 1, eps = 0, n = max(length(R0), length(ell)),
        init = c(1 - init.infected, init.infected),
        init.infected = .Machine[["double.neg.eps"]],
        weights = rep(c(1, 0), c(1L, n - 1L)),
        F = function (x) 1, Fargs = list(),
        H = identity, Hargs = list(),
        root = NULL, root.max = 1L, root.break = TRUE,
        aggregate = FALSE, skip.Y = FALSE, \dots)

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
\item{eps}{
  a non-negative number giving the the ratio of the mean time spent
  infectious and the mean life expectancy; \code{eps = 0} implies that
  life expectancy is infinite (that there are no deaths).}
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
\item{F, H}{
  functions returning a numeric vector of length 1 or of length equal
  that of the first formal argument.  The body must be symbolically
  differentiable with respect to the first formal argument; see
  \code{\link{D}}.}
\item{Fargs, Hargs}{
  lists of arguments passed to \code{F} or \code{H}.  In typical usage,
  \code{F} and \code{H} define parametric families of functions of one
  variable, and \code{Fargs} and \code{Hargs} supply parameter values.
  For example:
  \code{H = function(x, h) x^h}, \code{Hargs = list(h = 0.996)}.}
\item{root}{
  a function returning a numeric vector of length 1, with formal
  arguments \code{(tau, S, I, Y, dS, dI, dY, R0, ell)} (or a subset);
  otherwise, \code{NULL}.  Roots of this function in the interval from
  \code{from} to \code{to} are sought alongside the numerical solution.}
\item{root.max}{
  a positive integer giving a stopping condition for the root finder.
  Root finding continues until the count of roots found is
  \code{root.max}.}
\item{root.break}{
  a logical indicating if the solver should stop once \code{root.max}
  roots are found.  If \code{TRUE}, then the numerical solution ends
  at the last time point less than or equal to the last root.}
\item{aggregate}{
  a logical indicating if infected compartments should be aggregated.}
\item{skip.Y}{
  a logical indicating if solution of the equation for \eqn{Y} should
  not be attempted, e.g., because that equation seems stiffer than the
  rest.}
\item{\dots}{
  optional arguments passed to the solver, namely function
  \code{\link[deSolve]{lsoda}} in package \CRANpkg{deSolve}.}
\item{object}{
  an \R{} object inheriting from class \code{sir.aoi}, typically
  the value of a call to \code{sir.aoi}.}
\item{name}{
  a character string in \code{colnames(object)}.  Tail exponents of
  \eqn{V} are approximated, where, for example, \eqn{V = Y} if
  \code{name = "Y"} and \eqn{V = \sum_{j} I_{j}}{V = sum(I)}
  (prevalence) if \code{name = "I"}.}
\item{tol}{
  a positive number giving an upper bound on the relative change (from
  one time point to the next) in the slope of \eqn{\log(V)}{log(V)},
  defining time windows in which growth or decay of \eqn{V} is
  considered to be exponential.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{sir.aoi} and transitively from class \code{\link[=ts]{mts}},
storing the numerical solution.  Beneath the class, it is a
\code{T}-by-\code{(1+n+1)} numeric matrix of the form
\code{cbind(S, I, Y)}, \code{T <= length(seq(from, to, by))}.

If \code{root} is a function, then an attribute \code{root.info} of the
form \code{list(tau, state = cbind(S, I, Y))} stores the first \code{K}
roots of that function and the state of the system at each root,
\code{K <= root.max}.

If \code{aggregate = TRUE}, then infected compartments are aggregated
so that the number of columns named \code{I} is 1 rather than \code{n}.
This column stores prevalence, the proportion of the population that is
infected.  For convenience, there are 5 additional columns named
\code{I.E}, \code{I.I}, \code{foi}, \code{inc}, and \code{crv}.
These store the non-infectious and infectious components of prevalence
(so that \code{I.E + I.I = I}), the force of infection, incidence
(so that \code{foi * S = inc}), and the curvature of \code{Y}.

The method for \code{\link{summary}} returns a numeric vector of length
2 containing the left and right \dQuote{tail exponents} of the variable
\eqn{V} indicated by \code{name}, defined as the asymptotic values of
the slope of \eqn{\log(V)}{log(V)}.  \code{NaN} elements indicate that
a tail exponent cannot be approximated because the time window
represented by \code{object} does not cover enough of the tail, where
the meaning of \dQuote{enough} is set by \code{tol}.
}
\details{
The SIR equations with rates of transmission and recovery structured by
age of infection are
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} t
  &{} = \mu (1 - S) - (\textstyle\sum_{j} \beta_{j} F I_{j}) H(S) \\
  \text{d} & I_{    1} &{} / \text{d} t
  &{} = (\textstyle\sum_{j} \beta_{j} F I_{j}) H(S) - (\gamma_{1} + \mu) I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} t
  &{} = \gamma_{j} I_{j} - (\gamma_{j + 1} + \mu) I_{j + 1} \\
  \text{d} & R         &{} / \text{d} t
  &{} = \gamma_{n} I_{n} - \mu R
  \end{alignedat}
}{
  dS       /dt = mu * (1 - S) - sum(beta * F * I) * H(S)
  dI[    1]/dt = sum(beta * F * I) * H(S) - (gamma[1] + mu) * I[1]
  dI[j + 1]/dt = gamma[j] * I[j] - (gamma[j + 1] + mu) * I[j + 1]
  dR       /dt = gamma[n] * I[n] - mu * R
}
where
\eqn{S, I_{j}, R \ge 0}{min(S, I, R) >= 0},
\eqn{S + \sum_{j} I_{j} + R = 1}{sum(S, I, R) = 1},
\eqn{F} is a forcing function, and
\eqn{H} is a susceptible heterogeneity function.
In general, \eqn{F} and \eqn{H} are nonlinear.  In the standard
SIR equations, \eqn{F} is 1 and \eqn{H} is the identity function.

Nondimensionalization using parameters
\eqn{\mathcal{R}_{0,j} = \beta_{j}/(\gamma_{j} + \mu)}{R0 = beta/(gamma + mu)},
\eqn{\ell_{j} = (1/(\gamma_{j} + \mu))/t_{+}}{ell = (1/(gamma + mu))/tp}, and
\eqn{\varepsilon = t_{+}/(1/\mu)}{eps = tp/(1/mu)}
and independent variable
\eqn{\tau = t/t_{+}}{tau = t/tp},
where
\eqn{t_{+} = \sum_{j:\mathcal{R}_{0,j} > 0} 1/(\gamma_{j} + \mu)}{tp = sum(1/(gamma[R0 > 0] + mu))}
designates as a natural time unit the mean time spent infectious, gives
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} \tau
  &{} = \varepsilon (1 - S) - (\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) F I_{j}) H(S) \\
  \text{d} & I_{    1} &{} / \text{d} \tau
  &{} = (\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) F I_{j}) H(S) - (1/\ell_{1} + \varepsilon) I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} \tau
  &{} = (1/\ell_{j}) I_{j} - (1/\ell_{j+1} + \varepsilon) I_{j + 1} \\
  \text{d} & R         &{} / \text{d} \tau
  &{} = (1/\ell_{n}) I_{n} - \varepsilon R
  \end{alignedat}
}{
  dS       /dtau = eps * (1 - S) - sum((R0/ell) * F * I) * H(S)
  dI[    1]/dtau = sum((R0/ell) * F * I) * H(S) - (1/ell[1] + eps) * I[1]
  dI[j + 1]/dtau = (1/ell[j]) * I[j] - (1/ell[j + 1] + eps) * I[j + 1]
  dR       /dtau = (1/ell[n]) * I[n] - eps * R
}
\code{sir.aoi} works with the nondimensional equations, dropping the
last equation (which is redundant given
\eqn{R = 1 - S - \sum_{j} I_{j}}{R = 1 - sum(S, I)}) and augments the
resulting system of \eqn{1 + n} equations with a new equation
\deqn{
  \text{d}Y/\text{d}\tau = (\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) F I_{j}) (H(S) - 1/(\textstyle\sum_{j} \mathcal{R}_{0, j} F))
}{
  dY/tau = sum((R0/ell) * F * I) * (H(S) - 1/sum(R0 * F))
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

soln.peak <- sir.aoi(to = to, by = by, R0 = R0, ell = ell,
                     root = function (S, R0) sum(R0) * S - 1,
                     aggregate = TRUE)
str(soln.peak)

info.peak <- attr(soln.peak, "root.info")
to <- 4 * info.peak[["tau"]] # a more principled endpoint

soln <- sir.aoi(to = to, by = by, R0 = R0, ell = ell,
                aggregate = TRUE, atol = 0, rtol = 1e-12)
##                                ^^^^      ^^^^
## 'atol', 'rtol', ... are passed to deSolve::lsoda

head(soln)
tail(soln)

plot(soln) # dispatching stats:::plot.ts
plot(soln, log = "y")

(lamb <- summary(soln)) # left and right tail exponents

xoff <- function (x, k) x - x[k]
k <- c(16L, nrow(soln)) # c(1L, nrow(soln)) worse due to transient

plot(soln[, "Y"], log = "y", ylab = "Y")
abline(v = info.peak[["tau"]], h = info.peak[["state"]][, "Y"],
       lty = 2, lwd = 2, col = "red")
for (i in 1:2)
	lines(soln[k[i], "Y"] * exp(lamb[i] * xoff(time(soln), k[i])),
	      lty = 2, lwd = 2, col = "red")

wrap <-
function (root, ...)
	attr(sir.aoi(to = to, by = by, R0 = R0, ell = ell,
	             root = root, aggregate = TRUE, ...),
	     "root.info")[["tau"]]
Ymax <- info.peak[["state"]][, "Y"]

## NB: want *simple* roots, not *multiple* roots
L <- list(function (Y) (Y - Ymax * 0.5)  ,
          function (Y) (Y - Ymax * 0.5)^2,
          function (Y) (Y - Ymax      )  ,
          function (Y) (Y - Ymax      )^2)
lapply(L, wrap)

## NB: critical values can be attained more than once
L <- list(function (Y, dY)             Y - Ymax * 0.5,
          function (Y, dY) if (dY > 0) Y - Ymax * 0.5 else 1,
          function (Y, dY) if (dY < 0) Y - Ymax * 0.5 else 1)
lapply(L, wrap, root.max = 2L)

})
}
