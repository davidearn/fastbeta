\name{deconvolve}
\alias{deconvolve}
\title{Richardson-Lucy deconvolution}
\description{TODO}
%\details{}
\usage{
deconvolve(x, prob = 1, delay = 1,
           start, tol = 1, iter.max = 20L, complete = FALSE)
}
\arguments{
\item{x}{.}
\item{prob}{.}
\item{delay}{.}
\item{start}{.}
\item{tol}{.}
\item{iter.max}{.}
\item{complete}{.}
}
%\value{}
%\seealso{}
\references{
Goldstein, E., Dushoff, J., Ma, J., Plotkin, J. B., Earn, D. J. D.,
& Lipsitch, M. (2020).
Reconstructing influenza incidence by deconvolution of daily mortality
time series.
\emph{Proceedings of the National Academy of Sciences U. S. A.},
\emph{106}(51), 21825-21829.
\doi{10.1073/pnas.0902958106}
}
\examples{
beta <- function(t) 1e-05
nu   <- function(t) 0
mu   <- function(t) 0

S0 <- 5e+04
I0 <- 1e+03
constants <- c(gamma = 0.5, S0 = S0, I0 = I0, R0 = 1e+06 - S0 - I0)

set.seed(0)
n <- 250L
X <- sir(n, beta, nu, mu, constants, useCompiled = FALSE)
X # FIXME: many NA, apparently because adaptivetau ignores maxtau ... ??
}
