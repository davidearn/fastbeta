\name{fastbeta.matrix}
\alias{fastbeta.matrix}
\title{Calculate Coefficient Matrix for Iteration Step}
\description{
Calculates the coefficient matrix corresponding to one step of the
iteration carried out by \code{\link{fastbeta}}:
\preformatted{x <- c(1, E, I, R, S)
for (pos in seq_len(nrow(series) - 1L)) {
	L <- fastbeta.matrix(pos, series, ...)
	x <- L \%*\% x
}
}
}
\usage{
fastbeta.matrix(pos,
                series, sigma = gamma, gamma = 1, delta = 0,
                m = 1L, n = 1L)
}
\arguments{
\item{pos}{
  an integer indexing a row (but not the last row) of \code{series}.}
\item{series}{
  a \dQuote{multiple time series} object, inheriting from class
  \code{\link[=ts]{mts}}, with three columns storing (\dQuote{parallel},
  equally spaced) time series of incidence, births, and the per capita
  natural mortality rate, in that order.}
\item{sigma, gamma, delta}{
  non-negative numbers.  \code{m*sigma}, \code{n*gamma}, and
  \code{delta} are the rates of removal from each latent, infectious,
  and recovered compartment.}
\item{m}{
  a non-negative integer indicating a number of latent stages.}
\item{n}{
  a positive integer indicating a number of infectious stages.}
}
\value{
A lower triangular matrix of size \code{m + n + 3}.
}
%\details{}
%\seealso{}
%\references{}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
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

args <- c(list(pos = 1L, series = series),
          a[c("sigma", "gamma", "delta", "m", "n")])
str(args)

L <- do.call(fastbeta.matrix, args)
str(L)
symnum(L != 0)

})
}
