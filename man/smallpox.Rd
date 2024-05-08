\name{smallpox}
\alias{smallpox}
\docType{data}
\title{Smallpox Mortality in London, England, 1661-1930}
\description{
Time series of deaths due to smallpox, deaths due to all causes, and
births in London, England, from 1661 to 1930, as recorded in the
London Bills of Mortality and the Registrar General's Weekly Returns.
}
\usage{
data(smallpox, package = "fastbeta")
}
\format{
A data frame with 13923 observations of 5 variables:
\describe{
\item{from}{
  start date of the record.}
\item{nday}{
  length of the record, which is the number of days (typically 7) over
  which deaths and births were counted.}
\item{smallpox}{
  count of deaths due to smallpox.}
\item{allcauses}{
  count of deaths due to all causes.}
\item{births}{
  count of births.}
}
}
\source{
A precise description of the data set and its correspondence to the
original source documents is provided in the reference.

A script generating the \code{smallpox} data frame from a CSV file
accompanying the reference is available as
\code{\link{system.file}("scripts", "smallpox.R", package = "fastbeta")}.
}
\references{
Krylova, O. & Earn, D. J. D. (2020).
Patterns of smallpox mortality in London, England, over three centuries.
\emph{PLOS Biology},
\emph{18}(12), Article e3000506, 1-27.
\doi{10.1371/journal.pbio.3000506}
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
data(smallpox, package = "fastbeta")
str(smallpox)
table(smallpox[["nday"]]) # not all 7 days, hence:
plot(7 * smallpox / as.double(nday) ~ from, smallpox, type = "l")
}
