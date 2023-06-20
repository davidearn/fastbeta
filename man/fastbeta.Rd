\name{fastbeta}
\alias{fastbeta}
\title{Estimate a time-varying infectious disease transmission rate}
\description{
  Generates a discrete time estimate of the transmission rate
  \eqn{\beta(t)} from an equally spaced incidence time series
  and other data.
}
\details{
  \code{fastbeta} works by discretizing the SIR model
  \deqn{
    %% https://stat.ethz.ch/pipermail/r-devel/2023-June/082658.html
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
    X_{k} &\sim X(k) \qquad [X = S, I, \beta, \mu] \\
    Z_{k} &\sim \int_{k - 1}^{k} \beta(s) S(s) I(s) \, \text{d} s \\
    B_{k} &\sim \int_{k - 1}^{k} \nu(s) \, \text{d} s
    \end{aligned}
  }{
    X[k] ~ X(k)    [x = S, I, beta, mu]
    Z[k] ~ integrate(beta(s) S(s) I(s), k-1, k)
    B[k] ~ integrate(            nu(s), k-1, k)
  }
  and it is understood that indexing starts at \eqn{k = 0}
  (contrary to the \R{} convention).
}
\usage{
fastbeta(Z, B, mu, gamma, S0, I0)
}
\arguments{
  \item{Z}{a numeric vector specifying an incidence time series
    \eqn{(k, Z_{k})}{(k, Z[k])}.}
  \item{B}{a numeric vector specifying a births time series
    \eqn{(k, B_{k})}{(k, B[k])}.}
  \item{mu}{a numeric vector specifying a per capita natural mortality
    rate time series \eqn{(k, \mu_{k})}{(k, mu[k])}.}
  \item{gamma}{removal rate \eqn{\gamma}{gamma}.}
  \item{S0}{initial number of susceptible individuals \eqn{S_{0}}{S[0]}.}
  \item{I0}{initial number of infected individuals \eqn{I_{0}}{I[0]}.}
}
\value{
  A data frame with variables \code{S}, \code{I}, and \code{beta}
  corresponding elementwise to \code{Z}, \code{B}, and \code{mu}.
}
\seealso{\code{\link{ptpi}}.}
\references{
  Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
  Fast estimation of time-varying infectious disease transmission rates.
  \emph{PLoS Computational Biology},
  \emph{16}(9), Article e1008124, 1-39.
  \doi{10.1371/journal.pcbi.1008124}
}
%\examples{}
