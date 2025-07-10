\name{influenza}
\alias{influenza}
\title{Influenza and Pneumonia Mortality in Philadelphia, PA, 1918}
\description{
Time series of deaths due to influenza and pneumonia in Philadelphia, PA
from September 1, 1918 to December 31, 1918, as recorded in the
\dQuote{Special Tables of Mortality} of the U.S. Census Bureau.
}
\usage{
data(influenza, package = "fastbeta")
}
\format{
A data frame with 122 observations of 2 variables:
\describe{
\item{date}{
  date of the record.}
\item{deaths}{
  count of deaths due to influenza and pneumonia.}
}
}
\source{
Table 2 in the reference.
}
\references{
U.S. Census Bureau (1920).
\emph{Special Tables of Mortality from Influenza and Pneumonia:
  Indiana, Kansas, and Philadelphia, PA.}
U.S. Department of Commerce.
\url{https://www.census.gov/library/publications/1920/demo/1918-mortality-special-tables.html}
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
data(influenza, package = "fastbeta")
str(influenza)
plot(deaths ~ date, influenza)
}
\docType{data}
