\name{fastbeta}
\alias{fastbeta}
\title{Estimate a Time-Varying Infectious Disease Transmission Rate}
\description{
Generates a discrete approximation of a time-varying infectious disease
transmission rate from an equally spaced disease incidence time series
and other data.
}
\usage{
fastbeta(series, sigma = gamma, gamma = 1, delta = 0,
         init, m = length(init) - n - 2L, n = 1L, \dots)
}
\arguments{
\item{series}{
  a \dQuote{multiple time series} object, inheriting from class
  \code{\link[=ts]{mts}}, with three columns storing (\dQuote{parallel},
  equally spaced) time series of incidence, births, and the per capita
  natural mortality rate, in that order.}
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
\item{\dots}{
  optional arguments passed to \code{\link{deconvolve}}, if the first
  column of \code{series} represents \emph{observed} incidence rather
  than actual or estimated incidence.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, with \code{1+m+n+1+1} columns (named \code{S},
\code{E}, \code{I}, \code{R}, and \code{beta}) storing the result of the
iteration described in \sQuote{Details}.  It is completely parallel to
argument \code{series}, having the same \code{\link{tsp}} attribute.
}
\details{
The algorithm implemented by \code{fastbeta} is based on an SEIR model
with
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
It is derived by linearizing of the system of ordinary differential
equations
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
and substituting actual or estimated incidence and births for definite
integrals of \eqn{\lambda S} and \eqn{\nu}.  This procedure yields a
system of linear difference equations from which one recovers a discrete
approximation of \eqn{\beta}:
\deqn{
  \begin{alignedat}{17}
  &E_{t + 1}^{    1}
  &{} = {}& [(1 - \tfrac{1}{2} ( & m \sigma + \mu_{t})) & E_{t}^{    1} &       &              &            &           &       &               &   &{} + Z_{t + 1} &              ] &{} /
  [1 + \tfrac{1}{2} ( & m \sigma + \mu_{t + 1})] \\
  &E_{t + 1}^{i + 1}
  &{} = {}& [(1 - \tfrac{1}{2} ( & m \sigma + \mu_{t})) & E_{t}^{i + 1} &{} + {}& \tfrac{1}{2} & m \sigma ( & E_{t}^{i} &{} + {}& E_{t + 1}^{i} & ) &               &              ] &{} /
  [1 + \tfrac{1}{2} ( & m \sigma + \mu_{t + 1})] \\
  &I_{t + 1}^{    1}
  &{} = {}& [(1 - \tfrac{1}{2} ( & n \gamma + \mu_{t})) & I_{t}^{    1} &{} + {}& \tfrac{1}{2} & m \sigma ( & E_{t}^{m} &{} + {}& E_{t + 1}^{m} & ) &               &              ] &{} /
  [1 + \tfrac{1}{2} ( & n \gamma + \mu_{t + 1})] \\
  &I_{t + 1}^{j + 1}
  &{} = {}& [(1 - \tfrac{1}{2} ( & n \gamma + \mu_{t})) & I_{t}^{j + 1} &{} + {}& \tfrac{1}{2} & n \gamma ( & I_{t}^{j} &{} + {}& I_{t + 1}^{j} & ) &               &              ] &{} /
  [1 + \tfrac{1}{2} ( & n \gamma + \mu_{t + 1})] \\
  &R_{t + 1}
  &{} = {}& [(1 - \tfrac{1}{2} ( &   \delta + \mu_{t})) & R_{t}         &{} + {}& \tfrac{1}{2} & n \gamma ( & I_{t}^{n} &{} + {}& I_{t + 1}^{n} & ) &               &              ] &{} /
  [1 + \tfrac{1}{2} ( &   \delta + \mu_{t + 1})] \\
  &S_{t + 1}
  &{} = {}& [(1 - \tfrac{1}{2} ( &            \mu_{t})) & S_{t}         &{} + {}& \tfrac{1}{2} &   \delta ( & R_{t}     &{} + {}& R_{t + 1}     & ) &{} - Z_{t + 1} &{} + B_{t + 1}] &{} /
  [1 + \tfrac{1}{2} ( &            \mu_{t + 1})]
  \end{alignedat} \\
  \beta_{t} = (Z_{t} + Z_{t + 1}) / (2 S_{t} \sum_{j} I_{t}^{j})
}{
  E[t + 1,     1] = ( (1 - 0.5 * (m * sigma + mu[t])) * E[t, 1]                                             + Z[t + 1]            ) / ( 1 + 0.5 * (m * sigma + mu[t + 1]) )
  E[t + 1, i + 1] = ( (1 - 0.5 * (m * sigma + mu[t])) * E[t, i] + 0.5 * m * sigma * (E[t, i] + E[t + 1, i])                       ) / ( 1 + 0.5 * (m * sigma + mu[t + 1]) )
  I[t + 1,     1] = ( (1 - 0.5 * (n * gamma + mu[t])) * I[t, 1] + 0.5 * m * sigma * (E[t, m] + E[t + 1, m])                       ) / ( 1 + 0.5 * (n * gamma + mu[t + 1]) )
  I[t + 1, j + 1] = ( (1 - 0.5 * (n * gamma + mu[t])) * I[t, j] + 0.5 * n * gamma * (I[t, j] + I[t + 1, j])                       ) / ( 1 + 0.5 * (n * gamma + mu[t + 1]) )
  R[t + 1       ] = ( (1 - 0.5 * (    delta + mu[t])) * R[t   ] + 0.5 * n * gamma * (I[t, n] + I[t + 1, n])                       ) / ( 1 + 0.5 * (    delta + mu[t + 1]) )
  S[t + 1       ] = ( (1 - 0.5 * (            mu[t])) * S[t   ] + 0.5 *     delta * (R[t   ] + R[t + 1   ]) - Z[t + 1] + B[t + 1] ) / ( 1 + 0.5 * (            mu[t + 1]) )

  beta[t] = ( Z[t] + Z[t + 1] ) / ( 2 * S[t] * sum(I[t, ]) )
}
where we use the notation
\deqn{
  X_{t} \sim X(t) : X = S, E^{i}, I^{j}, R, Z, B, \mu, \beta \\
  \begin{aligned}
  Z(t) &= \int_{t - 1}^{t} \lambda(s) S(s) \, \text{d} s \\
  B(t) &= \int_{t - 1}^{t}          \nu(s) \, \text{d} s
  \end{aligned}
}{
  X[t   ] ~ X   (t) : X = S, R, Z, B, mu, beta
  X[t, .] ~ X[.](t) : X = E, I

  Z(t) = integrate(lambda * S, t - 1, t)
  B(t) = integrate(        nu, t - 1, t)
}
and it is understood that the independent variable \eqn{t} is a unitless
measure of time relative to the spacing of the substituted time series
of incidence and births.
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
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(   stats, pos = "package:base", verbose = FALSE)
library(graphics, pos = "package:base", verbose = FALSE)
library(   utils, pos = "package:base", verbose = FALSE)
}
if (requireNamespace("adaptivetau")) withAutoprint({

data(seir.ts02, package = "fastbeta")
a <- attributes(seir.ts02)
str(seir.ts02)
plot(seir.ts02)

## We suppose that we have perfect knowledge of incidence,
## births, and the data-generating parameters
series <- cbind(seir.ts02[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu") # FIXME: stats:::cbind.ts mangles dimnames

args <- c(list(series = series),
          a[c("sigma", "gamma", "delta", "init", "m", "n")])
str(args)

X <- do.call(fastbeta, args)
str(X)
plot(X)

plot(X[, "beta"], ylab = "Transmission rate")
lines(a[["beta"]](time(X)), col = "red") # the "truth"

})
}
