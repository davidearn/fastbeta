\name{sir.canonical}
\alias{sir.canonical}
\alias{summary.sir.canonical}
\title{Solve the Canonical SIR Equations}
\description{
Numerically integrates the canonical SIR equations with transmission
structured by age of infection.
}
\usage{
sir.canonical(from = 0, to = from + 1, by = 1,
              R0, ell, n = max(length(R0), length(ell)),
              init = c(1 - p, p), p = .Machine[["double.neg.eps"]],
              weights = rep(c(1, 0), c(1L, n - 1L)),
              root = c("none", "peak"), aggregate = FALSE, \dots)

\method{summary}{sir.canonical}(object, tol = 1e-6, \dots)
}
\arguments{
\item{from, to, by}{
  passed to \code{\link{seq.int}} in order to generate an increasing,
  equally spaced vector of time points in units of the mean time spent
  infectious.}
\item{R0}{
  a numeric vector of length \code{n} such that \code{sum(R0)} is the
  basic reproduction number and \code{R0[j]} is the contribution of
  infected compartment \code{j}.  \code{R0} of length 1 is recycled to
  length \code{n}.}
\item{ell}{
  a numeric vector of length \code{n} such that \code{ell[j]} is the
  ratio of the mean time spent in infected compartment \code{j} and the
  mean time spent infectious.  \code{ell} of length 1 is recycled to
  length \code{n}.  \code{ell/sum(ell[R0 > 0])} is used internally.}
\item{n}{
  a positive integer giving the number of infected compartments.}
\item{init}{
  a numeric vector of length 2 giving initial susceptible and infected
  proportions.}
\item{p}{
  a number in \eqn{(0, 1]} used only to define \code{init} when
  \code{init} is unset.}
\item{weights}{
  a numeric vector of length \code{n} containing non-negative weights,
  defining the initial distribution of infected individuals among the
  infected compartments.  By default, all infected individuals occupy
  the first compartment.}
\item{root}{
  a character string determining what is returned.  \code{"none"}: the
  numerical solution as a multiple time series object.  \code{"peak"}:
  information about the time and state when \eqn{Y} attains a local
  maximum.}
\item{aggregate}{
  a logical indicating if infected compartments should be aggregated
  when \code{root = "none"}.}
\item{\dots}{
  optional arguments passed to \code{\link[deSolve]{lsoda}}.}
\item{object}{
  an \R{} object inheriting from class \code{sir.canonical}, typically
  the value of a call to \code{sir.canonical}.}
\item{tol}{
  a positive number giving an upper bound on the relative change (from
  one time point to the next) in the slope of log prevalence, defining
  time windows in which growth or decay of prevalence is considered to
  be exponential.}
}
\value{
If \code{root = "none"},
a \dQuote{multiple time series} object, inheriting from class
\code{sir.canonical} and transitively from class \code{\link[=ts]{mts}}.
Beneath the class, it is a
\code{length(seq(from, to, by))}-by-\code{(1+n+1)} numeric matrix with
columns \code{S}, \code{I}, and \code{Y}. % FIXME: aggregate=TRUE

If \code{root = "peak"},
a numeric vector of length 6 of the form
\code{c(tau, S, I, Y, I.E, I.I)}
containing the time in units of the mean time spent infectious at which
\eqn{Y} attains a local maximum and the values at that time of
\eqn{S},
\eqn{\sum_{j} I_{j}}{sum(I)},
\eqn{Y},
\eqn{\sum_{j:\mathcal{R}_{0,j} = 0} I_{j}}{sum(I[R0 == 0])}, and
\eqn{\sum_{j:\mathcal{R}_{0,j} > 0} I_{j}}{sum(I[R0 > 0])}.
Attributes \code{I} and \code{curvature} store the corresponding values
of \eqn{I_{j}}{I[j]} and \eqn{Y''}.  If \eqn{Y} does not attain a local
maximum between times \code{from} and \code{to}, then the value is
\code{NULL}.

The method for \code{\link{summary}} returns a numeric vector of length
2 containing the left and right \dQuote{tail exponents}, defined as the
asymptotic values of the slope of log prevalence.  \code{NaN} elements
indicate that a tail exponent cannot be approximated from the prevalence
time series represented by \code{object}, because the time window does
not cover enough of the tail, where the meaning of \dQuote{enough} is
set by \code{tol}.
}
\details{
The standard SIR equations with age-structured transmission are
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} t
  &{} = -\textstyle\sum_{j} (\beta_{j}/N) S I_{j} \\
  \text{d} & I_{    1} &{} / \text{d} t
  &{} = \textstyle\sum_{j} (\beta_{j}/N) S I_{j} - \gamma_{1} I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} t
  &{} = \gamma_{j} I_{j} - \gamma_{j + 1} I_{j + 1} \\
  \text{d} & R         &{} / \text{d} t
  &{} = \gamma_{n} I_{n}
  \end{alignedat}
}{
  dS       /dt = -sum((beta/N) * S * I)
  dI[    1]/dt = sum((beta/N) * S * I) - gamma[1] * I[1]
  dI[j + 1]/dt = gamma[j] * I[j] - gamma[j + 1] * I[j + 1]
  dR       /dt = gamma[n] * I[n]
}
where \eqn{N = S + \sum_{j} I_{j} + R}{N = sum(S, I, R)} is the
(constant, positive) population size.  Without loss of generality,
the canonical SIR equations with age-structured transmission are
obtained by substituting \eqn{N = 1},
\eqn{\mathcal{R}_{0,j} = \beta_{j}/\gamma_{j}}{R0 = beta/gamma}, and
\eqn{\ell_{j} = (1/\gamma_{j})/\sum_{j:\mathcal{R}_{0,j} > 0} (1/\gamma_{j})}{ell = (1/gamma)/sum(1/gamma[R0 > 0])}
and choosing
\eqn{\tau = t/\sum_{j:\mathcal{R}_{0,j} > 0} (1/\gamma_{j})}{tau = t/sum(1/gamma[R0 > 0])}
as the independent variable:
\deqn{
  \begin{alignedat}{4}
  \text{d} & S         &{} / \text{d} \tau
  &{} = -\textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) S I_{j} \\
  \text{d} & I_{    1} &{} / \text{d} \tau
  &{} = \textstyle\sum_{j} (\mathcal{R}_{0,j}/\ell_{j}) S I_{j} - (1/\ell_{1}) I_{1} \\
  \text{d} & I_{j + 1} &{} / \text{d} \tau
  &{} = (1/\ell_{j}) I_{j} - (1/\ell_{j+1}) I_{j + 1} \\
  \text{d} & R         &{} / \text{d} \tau
  &{} = (1/\ell_{n}) I_{n} \\
  \end{alignedat}
}{
  dS       /dtau = -sum((R0/ell) * S * I)
  dI[    1]/dtau = sum((R0/ell) * S * I) - (1/ell[1]) * I[1]
  dI[j + 1]/dtau = (1/ell[j]) * I[j] - (1/ell[j + 1]) * I[j + 1]
  dR       /dtau = (1/ell[n]) * I[n]
}
\code{sir.canonical} drops the last equation, which is redundant as
\eqn{R = 1 - S - \sum_{j} I_{j}}{R = 1 - sum(S, I)}, and augments the
resulting system of \eqn{1 + n} equations with a new equation
\deqn{
  \text{d}Y/\text{d}\tau = (\sum_{j} \mathcal{R}_{0, j} S - 1) \sum_{j:\mathcal{R}_{0,j} > 0} I_{j}
}{
  dY/tau = (sum(R0) * S - 1) * sum(I[R0 > 0])
}
due to the usefulness of the solution \eqn{Y} in applications.
}
\note{
\code{sir.canonical} is not a special case of \code{\link{sir}} nor a
generalization.  The two functions were developed independently and for
different purposes: \code{sir.canonical} to validate analytical results
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

to <- 100; by <- 0.01; R0 <- 16; ell <- 1/3

peak <- sir.canonical(to = to, by = by, R0 = R0, ell = ell,
                      root = "peak")
peak

to <- 4 * peak[["tau"]] # a more principled endpoint

soln <- sir.canonical(to = to, by = by, R0 = R0, ell = ell,
                      aggregate = TRUE)
head(soln)

plot(soln) # dispatching stats:::plot.ts

plot(soln[, "Y"])
abline(v = peak[["tau"]], h = peak[["Y"]],
       lty = 2, lwd = 2, col = "red")

xoff <- function (x, k) x - x[k]
prev <- soln[, "I"]
lamb <- summary(soln)
k <- c(16L, nrow(soln)) # c(1L, nrow(soln)) suffers due to transient
plot(prev, log = "y", ylab = "Prevalence")
for (i in 1:2)
	lines(prev[k[i]] * exp(lamb[i] * xoff(time(prev), k[i])),
	      lty = 2, lwd = 2, col = "red")

})
}
