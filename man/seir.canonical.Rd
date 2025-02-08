\name{seir.canonical}
\alias{seir.canonical}
\title{Solve the Canonical SEIR Equations}
\description{
Numerically integrates the canonical SEIR equations, a special case of
the more general SEIR equations handled by \code{\link{seir}}; see
\sQuote{Details}.
}
\usage{
seir.canonical(from = 0, to = from + 1, by = 1,
               R0, ell = if (m == 0L) 0 else (2 * n)/(3 * n + 1),
               m = 1L, n = 1L, init = c(1 - p, p), p = 0x1p-64,
               weights = rep(c(1, 0), c(1L, m + n - 1L)),
               root = c("none", "peak"), \dots)
}
\arguments{
\item{from, to, by}{
  passed to \code{\link{seq.int}} in order to generate an increasing,
  equally spaced vector of time points in units of the generation
  interval.}
\item{R0}{
  a positive number; the basic reproduction number.}
\item{ell}{
  a number in \eqn{[0,1)}, equal to zero if and only if \code{m=0}; the
  ratio of the mean latent period and mean generation interval.}
\item{m}{
  a non-negative integer indicating a number of latent stages.}
\item{n}{
  a positive integer indicating a number of infectious stages.}
\item{init}{
  a numeric vector of length 2 giving initial susceptible and infected
  proportions.}
\item{p}{
  a number in \eqn{(0,1]} defining the default value of \code{init}.  It
  is an error to set both of \code{init} and \code{p}.}
\item{weights}{
  a numeric vector of length \code{m+n} containing non-negative weights,
  defining the initial distribution of infected individuals among the
  latent and infectious stages.  By default, all infected individuals
  occupy the first stage.}
\item{root}{
  a character string.  \code{"none"}: Return the numerical solution as a
  \dQuote{multiple time series} object.  \code{"peak"}: Return the time
  in units of the generation interval at which \eqn{Y} attains a local
  maximum, \code{NaN} if \eqn{Y'} has no root between \code{from} and
  \code{to}.}
\item{\dots}{
  optional arguments passed to \code{\link[deSolve]{lsoda}}.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}.  Beneath the class, it is a
\code{length(seq(from, to, by))}-by-\code{(1+m+n+1)} numeric matrix
with columns \code{S}, \code{E}, \code{I}, and \code{Y}.
}
\details{
The general SEIR equations handled by \code{\link{seir}} are
\deqn{
  \begin{alignedat}{5}
  \text{d} & S         &{} / \text{d} t
  &{} = {}& \nu(t) + \delta R - (\beta(t) \sum_{j} I^{j} + \mu(t)) S \\
  \text{d} & E^{    1} &{} / \text{d} t
  &{} = {}& \beta(t) S \sum_{j} I^{j} - (m \sigma + \mu(t)) E^{1} \\
  \text{d} & E^{i + 1} &{} / \text{d} t
  &{} = {}& m \sigma E^{i} - (m \sigma + \mu(t)) E^{i + 1} \\
  \text{d} & I^{    1} &{} / \text{d} t
  &{} = {}& m \sigma E^{m} - (n \gamma + \mu(t)) I^{1} \\
  \text{d} & I^{j + 1} &{} / \text{d} t
  &{} = {}& n \gamma I^{j} - (n \gamma + \mu(t)) I^{j + 1} \\
  \text{d} & R         &{} / \text{d} t
  &{} = {}& n \gamma I^{n} - (\delta + \mu(t)) R
  \end{alignedat}
}{
  dS       /dt = nu(t) + delta * R - (beta(t) * sum(I) + mu(t)) * S
  dE[    1]/dt = beta(t) * S * sum(I) - (m * sigma + mu(t)) * E[1]
  dE[i + 1]/dt = m * sigma * E[i] - (m * sigma + mu(t)) * E[i + 1]
  dI[    1]/dt = m * sigma * E[m] - (n * gamma + mu(t)) * I[1]
  dI[j + 1]/dt = n * gamma * I[j] - (n * gamma + mu(t)) * I[j + 1]
  dR       /dt = n * gamma * I[n] - (delta + mu(t)) * R
}
The canonical SEIR equations are obtained by substituting
\eqn{\beta(t) = \mathcal{R}_{0} \gamma}{beta(t) = R0 * gamma},
\eqn{\nu(t) = \mu(t) = \delta = 0},
\eqn{1/\sigma = \ell g}{1/sigma = ell * g}, and
\eqn{1/\gamma = (1 - \ell) g/h}{1/gamma = (1 - ell) * g/h}, where
\eqn{h = (n + 1)/(2 n)}{h = (n + 1)/(2 * n)}, then choosing
\eqn{\tau = t/g} as an independent variable:
\deqn{
  \begin{alignedat}{5}
  \text{d} & S         &{} / \text{d} \tau
  &{} = {}& -(h \mathcal{R}_{0}/(1 - \ell)) S \sum_{j} I^{j} \\
  \text{d} & E^{    1} &{} / \text{d} \tau
  &{} = {}& (h \mathcal{R}_{0}/(1 - \ell)) S \sum_{j} I^{j} - (m/\ell) E^{1} \\
  \text{d} & E^{i + 1} &{} / \text{d} \tau
  &{} = {}& (m/\ell) E^{i} - (m/\ell) E^{i + 1} \\
  \text{d} & I^{    1} &{} / \text{d} \tau
  &{} = {}& (m/\ell) E^{m} - (h n/(1 - \ell)) I^{1} \\
  \text{d} & I^{j + 1} &{} / \text{d} \tau
  &{} = {}& (h n/(1 - \ell)) I^{j} - (h n/(1 - \ell)) I^{j + 1} \\
  \text{d} & R         &{} / \text{d} \tau
  &{} = {}& (h n/(1 - \ell)) I^{n}
  \end{alignedat}
}{
  dS       /dtau = -(h * R0/(1 - ell)) * S * sum(I)
  dE[    1]/dtau = (h * R0/(1 - ell)) * S * sum(I) - (m/ell) * E[1]
  dE[i + 1]/dtau = (m/ell) * E[i] - (m/ell) * E[i + 1]
  dI[    1]/dtau = (m/ell) * E[m] - (h * n/(1 - ell)) * I[1]
  dI[j + 1]/dtau = (h * n/(1 - ell)) * I[j] - (h * n/(1 - ell)) * I[j + 1]
  dR       /dtau = (h * n/(1 - ell)) * I[n]
}
The constraint \eqn{\mu(t) = \nu(t)} implies that the population size
\eqn{N = S + \sum_{i} E^{i} + \sum_{j} I^{j} + R} is constant, hence,
without loss of generality, \code{seir.canonical} assumes \eqn{N = 1}
and drops the last equation.  The resulting system of \eqn{1 + m + n}
equations is augmented with a final equation
\deqn{
  \text{d}Y/\text{d}\tau = ((\mathcal{R}_{0}/(1 - \ell)) S - 1) \sum_{j} I^{j}
}{
  dY/tau = ((R0/(1 - ell)) * S - 1) * sum(I)
}
due to the usefulness of the solution \eqn{Y} in applications.
}
\examples{
## TODO
}
