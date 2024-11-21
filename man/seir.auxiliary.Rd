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
Currently, none of the functions documented here are vectorized.
Arguments must have length 1.
}
\seealso{\code{\link{seir}}.}
%\examples{}
