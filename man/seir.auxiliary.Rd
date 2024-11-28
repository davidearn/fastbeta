\name{seir.auxiliary}
\alias{seir.auxiliary}
\alias{seir.R0}
\alias{seir.ee}
\alias{seir.jacobian}
\title{Auxiliary Functions for the SEIR Model without Forcing}
\description{
Calculate the basic reproduction number, endemic equilibrium,
and Jacobian matrix of the SEIR model without forcing.
}
\usage{
seir.R0      (beta, nu = 0, mu = 0, sigma = 1, gamma = 1, delta = 0,
              m = 1L, n = 1L, N = 1)
seir.ee      (beta, nu = 0, mu = 0, sigma = 1, gamma = 1, delta = 0,
              m = 1L, n = 1L, N = 1)
seir.jacobian(beta, nu = 0, mu = 0, sigma = 1, gamma = 1, delta = 0,
              m = 1L, n = 1L)
}
\arguments{
\item{beta, nu, mu, sigma, gamma, delta}{
  non-negative numbers.  \code{beta}, \code{nu}, and \code{mu} are the
  rates of transmission, birth, and natural death.  \code{m*sigma},
  \code{n*gamma}, and \code{delta} are the rates of removal from each
  latent, infectious, and recovered compartment.}
\item{m}{
  a non-negative integer indicating a number of latent stages.}
\item{n}{
  a positive integer indicating a number of infectious stages.}
\item{N}{
  a non-negative number indicating a population size for the
  \code{(nu == 0 && mu == 0)} case.}
}
\value{
\code{seir.R0} returns a numeric vector of length 1.  \code{seir.ee}
returns a numeric vector of length \code{1+m+n+1}.  \code{seir.jacobian}
returns a function of one argument \code{x} (which must be a numeric
vector of length \code{1+m+n+1}) whose return value is a square numeric
matrix of size \code{length(x)}.
}
\details{
If \eqn{\mu, \nu = 0}{mu = 0 and nu = 0}, then the basic reproduction
number is computed as
\deqn{
  \mathcal{R}_{0} = N \beta / \gamma
}{
  R_0 = N * beta / gamma
}
and the endemic equilibrium is computed as
\deqn{
  \begin{bmatrix}
  S^{\hphantom{1}} \\
  E^{i} \\
  I^{j} \\
  R^{\hphantom{1}}
  \end{bmatrix}
  =
  \begin{bmatrix}
  \gamma / \beta \\
  w \delta / (m \sigma) \\
  w \delta / (n \gamma) \\
  w
  \end{bmatrix}
}{
  S    = gamma / beta
  E[i] = w * delta / (m * sigma)
  I[j] = w * delta / (m * sigma)
  R    = w
}
where \eqn{w} is chosen so that the sum is \eqn{N}.

If \eqn{\mu, \nu > 0}{mu > 0 and nu > 0}, then the basic reproduction
number is computed as
\deqn{
  \mathcal{R}_{0} = \nu \beta a^{-m} (1 - b^{-n}) / \mu^{2}
}{
  R_0 = nu * beta * a^-m * (1 - b^-n) / mu^2
}
and the endemic equilibrium is computed as
\deqn{
  \begin{bmatrix}
  S^{\hphantom{1}} \\
  E^{i} \\
  I^{j} \\
  R^{\hphantom{1}}
  \end{bmatrix}
  =
  \begin{bmatrix}
  \mu a^{m} / (\beta (1 - b^{-n})) \\
  w a^{m - i} b^{n} (\delta + \mu) / (m \sigma) \\
  w b^{n - j} (\delta + \mu) / (n \gamma) \\
  w
  \end{bmatrix}
}{
  S    = mu * a^m / (beta * (1 - b^-n))
  E[i] = w * a^(m - i) * b^n * (delta + mu) / (m * sigma)
  I[j] = w * b^(n - j)       * (delta + mu) / (n * gamma)
  R    = w
}
where \eqn{w} is chosen so that the sum is \eqn{\nu / \mu}{nu / mu},
the population size at equilibrium, and
\eqn{a = 1 + \mu / (m \sigma)}{a = 1 + mu / (m * sigma)} and
\eqn{b = 1 + \mu / (n \gamma)}{b = 1 + mu / (n * gamma)}.

Currently, none of the functions documented here are vectorized.
Arguments must have length 1.
}
\seealso{
\code{\link{seir}}, for the system of ordinary differential equations
on which these computations are predicated.
}
%\examples{}
