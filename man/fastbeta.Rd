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
and computing the recursion
\deqn{
  \begin{aligned}
  S_{k} &= \frac{(1 - \frac{1}{2} \mu_{k - 1}) S_{k - 1} + B_{k} - Z_{k}}{
    1 + \frac{1}{2} \mu_{k}} \\
  I_{k} &= \frac{(1 - \frac{1}{2} (\gamma + \mu_{k - 1})) I_{k - 1} + Z_{k}}{
    1 + \frac{1}{2} (\gamma + \mu_{k})} \\
  \beta_{k} &= \frac{Z_{k} + Z_{k + 1}}{2 S_{k} I_{k}}
  \end{aligned}
}{
  S[k] = ((1 - 0.5 * mu[k-1]) * S[k-1] + B[k] - Z[k]) / (1 + 0.5 * mu[k])
  I[k] = ((1 - 0.5 * (gamma + mu[k-1])) * I[k-1] + Z[k]) / (1 + 0.5 * (gamma + mu[k]))
  beta[k] = (Z[k] + Z[k+1]) / (2 * S[k] * I[k])
}
where
\deqn{
  \begin{aligned}
  X_{k} &\sim X(k) \quad [X = S, I, \beta, \mu] \\
  Z_{k} &\sim \int_{k - 1}^{k} \beta(s) S(s) I(s) \, \text{d} s \\
  B_{k} &\sim \int_{k - 1}^{k} \nu(s) \, \text{d} s
  \end{aligned}
}{
  X[k] ~ X(k)    [X = S, I, beta, mu]
  Z[k] ~ integrate(beta(s) S(s) I(s), k-1, k)
  B[k] ~ integrate(            nu(s), k-1, k)
}
and it is understood that indexing starts at \eqn{k = 0}
(contrary to the \R{} convention).
}
\usage{
fastbeta(data = ts(cbind(Z, B, mu), start = 0),
         Z, B, mu, gamma, S0, I0)
}
\arguments{
\item{data}{a \dQuote{multiple time series} object, inheriting
  from class \code{\link[=ts]{mts}}.  The first three columns
  should define (\dQuote{parallel}, equally spaced) time series
  of incidence, births, and the per capita natural mortality
  rate, respectively.}
\item{Z, B, mu}{numeric vectors used to construct \code{data}
  if it is missing and unused otherwise.}
\item{gamma}{removal rate.}
\item{S0}{the initial number of susceptible individuals.}
\item{I0}{the initial number of infected individuals.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, parallel to \code{data}.  It has three
columns, \code{S}, \code{I}, and \code{beta}; these store the
result of the recursion described in \sQuote{Details}.
}
%\seealso{}
\references{
Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
Fast estimation of time-varying infectious disease transmission rates.
\emph{PLoS Computational Biology},
\emph{16}(9), Article e1008124, 1-39.
\doi{10.1371/journal.pcbi.1008124}
}
\examples{
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

zzz <- fastbeta(data = ts(cbind(X[, "Z"], nu(0), mu(0)), start = 0),
                gamma = par[["gamma"]],
                S0 = par[["S0"]],
                I0 = par[["I0"]])
}
