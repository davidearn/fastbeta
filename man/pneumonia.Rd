\name{pneumonia}
\alias{pneumonia}
\title{Pneumonia and Influenza Mortality in Philadelphia, PA, 1918}
\description{
Time series of deaths due to pneumonia and influenza in Philadelphia, PA
from September 1, 1918 to December 31, 1918, as recorded in the
\dQuote{Special Tables of Mortality} of the U.S. Census Bureau.
}
\usage{
data(pneumonia, package = "fastbeta")
}
\format{
A named list with 2 components, \code{series} and \code{delay}.
\code{series} is a data frame with 122 rows and 2 variables:
\describe{
\item{date}{
  date of the record.}
\item{deaths}{
  count of deaths due to influenza and pneumonia.}
}
\code{delay} is a data frame with 64 rows and 3 variables:
\describe{
\item{nday}{
  number of days from infection to death.}
\item{goldstein, gpg}{
  probabilities, not summing to 1 due to rounding and truncation; see
  \sQuote{Source}.}
}
}
\source{
A script generating the \code{pneumonia} object is available as
\code{\link{system.file}("scripts", "pneumonia.R", package = "fastbeta")}.

\code{series} is obtained from Table 2 in the first reference.

\code{delay} is obtained from the remaining references.  Component
\code{goldstein} is obtained from Figure 1 in the Supporting Information
of Goldstein et al. (2009).  Component \code{gpg} is obtained from the
convolution of two gamma distributions, one for the time from infection
to symptom onset fitted to Figure 1 in Moser et al. (1979) and another
for the time from symptom onset to death fitted to Chart 2 in Keeton &
Cushman (1918).
}
\references{
U.S. Census Bureau (1920).
\emph{Special Tables of Mortality from Influenza and Pneumonia:
  Indiana, Kansas, and Philadelphia, PA.}
U.S. Department of Commerce.
\url{https://www.census.gov/library/publications/1920/demo/1918-mortality-special-tables.html}

Goldstein, E., Dushoff, J., Ma, J., Plotkin, J. B., Earn, D. J. D., &
Lipsitch, M. (2009).
Reconstructing influenza incidence by deconvolution of daily mortality
time series.
\emph{Proceedings of the National Academy of Sciences U. S. A.},
\emph{106}(51), 21825-21829.
\doi{10.1073/pnas.0902958106}

Moser, M. R., Bender, T. R., Margolis, H. S., Noble, G. R.,
Kendal, A. P., & Ritter, D. G. (1979).
An outbreak of influenza aboard a commercial airliner.
\emph{Americal Journal of Epidemiology},
\emph{110}(1), 1-6.
\doi{10.1093/oxfordjournals.aje.a112781}

Keeton, R. W. & Cushman, A. B. (1918).
The influenza epidemic in Chicago: the disease as a type of toxemic
shock.
\emph{Journal of the Americal Medical Association}.
\emph{71}(24), 1962-1967.
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(graphics, pos = "package:base", verbose = FALSE)
library(   utils, pos = "package:base", verbose = FALSE)
}
data(pneumonia, package = "fastbeta")
str(pneumonia)

plot(deaths ~ date, pneumonia$series, xlab = "1918")

plot(goldstein/sum(goldstein) ~ nday, pneumonia$delay, type = "o",
     lty = 2, pch = 1, xlab = "days", ylab = "probability")
lines(gpg/sum(gpg) ~ nday, pneumonia$delay, type = "o",
      lty = 1, pch = 16)
}
\docType{data}
