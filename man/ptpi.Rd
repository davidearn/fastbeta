\name{ptpi}
\alias{ptpi}
\title{Peak to peak iteration}
\description{
  TODO
}
\usage{
ptpi(cases, births, mu,
     start, a = 1L, b = length(cases),
     tol = 1e-06, iter.max = 20L, complete = FALSE)
}
\arguments{
  \item{cases}{a numeric vector specifying an incidence time series.}
  \item{births}{a numeric vector specifying a births time series.}
  \item{mu}{a numeric vector specifying a per capita natural mortality
    rate time series.}
  \item{start}{starting value for the number of susceptible individuals
    at index \code{a}.}
  \item{a}{0-based index of the first peak in the incidence time series.}
  \item{b}{0-based index of the last peak in the incidence time series
    that is in phase with the first.}
  \item{tol}{a tolerance indicating a stopping condition.  }
  \item{iter.max}{maximum number of iterations.}
  \item{complete}{a logical flag indicating if the result should
    preserve the number of susceptible individuals from index \code{a}
    to index \code{b} in each iteration.}
}
%\details{}
%\value{}
%\seealso{}
%\references{}
%\examples{}
