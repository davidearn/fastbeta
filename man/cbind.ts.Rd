\name{cbind.ts}
\alias{cbind.ts}
\title{Combine Time Series Objects}
\description{
A replacement for the S3 method registered by package \pkg{stats} for
generic function \code{\link{cbind}} and class \code{\link{ts}}.  It
sets column names following the rules employed by the internal default
method for \code{cbind}.  It exists to allow users to work around
\PR{18583}, which shows that the method in package \pkg{stats} employs
different and rather inconsistent rules.  This function must be called
directly, as it is not registered as a method for \code{cbind}.
}
\usage{
\method{cbind}{ts}(\dots, deparse.level = 1)
}
\arguments{
\item{\dots}{
  vectors (including matrices), at least one inheriting from class
  \code{\link{ts}}.}
\item{deparse.level}{
  an integer (0, 1, or 2) controlling how column names are determined
  for untagged arguments that are not matrices, following the internal
  default method for \code{\link{cbind}}.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}.
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
n <- 3L
x <- matrix(0, n, n, dimnames = list(NULL, LETTERS[seq_len(n)]))
y <- seq_len(n)
tsx <- ts(x)
tsy <- ts(y)
`~` <- identity
for (k in 0L:2L) {
	cat(sprintf("k = \%d:\n\n\n", k))
	withAutoprint({
	try(colnames(cbind   (   x,    y, deparse.level = k)))
	try(colnames(cbind   ( tsx,  tsy, deparse.level = k)))
	try(colnames(cbind.ts( tsx,  tsy, deparse.level = k)))
	try(colnames(cbind   (~  x, ~  y, deparse.level = k)))
	try(colnames(cbind   (~tsx, ~tsy, deparse.level = k)))
	try(colnames(cbind.ts(~tsx, ~tsy, deparse.level = k)))
	})
	cat("\n\n")
}
rm(`~`)
}
