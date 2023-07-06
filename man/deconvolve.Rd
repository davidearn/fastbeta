\name{deconvolve}
\alias{deconvolve}
\title{Richardson-Lucy deconvolution}
\description{TODO}
%\details{}
\usage{
deconvolve(x, prob = 1, delay = 1,
           start = 0.125 + c(double(length(delay) - 1L), x),
           tol = 1, iter.max = 20L, complete = FALSE)
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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
h <- function(x, a = 1, b = 1, c = 0) a * exp(-b * (x - c)^2)

set.seed(0)
n <- 250L
d <- 100L
delay <- diff(pgamma(0L:(d + 1L), 50))
x <- rpois(n, h(seq(0, 30, length.out = n), 100, 0.1, 10))
y <- tabulate(rep.int(1L:n, x) +
              sample(0L:d, size = sum(x), replace = TRUE, prob = delay),
              n)
z <- deconvolve(as.double(y), delay = delay)$dx[-(1L:d)]
matplot(cbind(x, y, z))
}
