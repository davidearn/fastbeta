\name{ptpi}
\alias{ptpi}
\title{Peak to peak iteration}
\description{
TODO
}
\usage{
ptpi(series, start, a = 0L, b = nrow(series) - 1L,
     tol = 1e-06, iter.max = 20L, complete = FALSE)
}
\arguments{
\item{series}{a \dQuote{multiple time series} object, inheriting from
  class \code{\link[=ts]{mts}}, with three columns storing
  (\dQuote{parallel}, equally spaced) time series of incidence, births,
  and the per capita natural mortality rate, in that order.}
\item{start}{starting value for the number of susceptible individuals at
  index \code{a}.}
\item{a}{0-index of the first peak in the incidence time series.}
\item{b}{0-index of the last peak in the incidence time series that is
  in phase with the first.}
\item{tol}{a tolerance indicating a stopping condition, namely
  \eqn{\left|\frac{S_{b} - S_{a}}{S_{a}}\right| < \mathtt{tol}}{
    abs((S[b] - S[a]) / S[a]) < tol}.}
\item{iter.max}{maximum number of iterations.}
\item{complete}{a logical flag indicating if the result should preserve
  the number of susceptible individuals from index \code{a} to index
  \code{b} in each iteration.}
}
%\details{}
\value{
A list with elements:
\item{value}{the value of \eqn{S_{0}}{S[0]} back-calculated from
  the value of \eqn{S_{b}}{S[b]} obtained in the last iteration.}
\item{delta}{the relative difference between the last obtained
  values of \eqn{S_{a}}{S[a]} and \eqn{S_{b}}{S[b]}.}
\item{iter}{the number of iterations performed.}
\item{X}{if \code{complete = TRUE}, then a matrix with
  \code{b-a+1} rows and \code{iter} columns storing the
  value of \eqn{S_{k}}{S[k]} obtained in each iteration,
  for \eqn{k = a,\ldots,b}.  Otherwise, \code{NULL}.}
}
%\seealso{}
\references{
Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
Fast estimation of time-varying infectious disease transmission rates.
\emph{PLoS Computational Biology},
\emph{16}(9), Article e1008124, 1-39.
\doi{10.1371/journal.pcbi.1008124}
}
%\examples{}
