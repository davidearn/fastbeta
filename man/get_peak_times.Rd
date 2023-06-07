\name{get_peak_times}
\alias{get_peak_times}
\title{Find times of peaks in periodic time series}
\usage{
get_peak_times(x, period, bw_mavg, bw_peakid)
}
\arguments{
\item{x}{Numeric vector. A periodic, equally spaced time series.
Time points are taken to be \code{seq_along(x)}.}

\item{period}{Numeric scalar. The period of \code{x} in units of the
observation interval.}

\item{bw_mavg}{Numeric scalar. Bandwidth for the central moving
average applied to \code{x}. \code{x_mavg[i]} will be the mean of the set
of observations \code{x[j]} with
\ifelse{latex}{\out{$|i - j| \leq bandwidth$}}{\ifelse{html}{\out{&vert;<i>i</i> &minus; <i>j</i>&vert; &le; bandwidth}}{|i - j| <= bandwidth}}.}

\item{bw_peakid}{Numeric scalar. Bandwidth for peak identification.
\code{i} will be a peak time if and only if \code{x_mavg[i] > x_mavg[j]}
for all \code{j} with
\ifelse{latex}{\out{$0 < |i - j| \leq bandwidth$}}{\ifelse{html}{\out{0 &lt; &vert;<i>i</i> &minus; <i>j</i>&vert; &le; bandwidth}}{0 < |i - j| <= bandwidth}}.}
}
\value{
A list containing:

\describe{
\item{\code{x_mavg}}{Numeric vector. The result of applying a
\code{2 * bw_mavg + 1} point central moving average to \code{x}.
}
\item{\code{all}}{Numeric vector. A subset of \code{seq_along(x)}
listing the times of all identified peaks.
}
\item{\code{phase}}{Numeric vector. A subset of \code{seq_along(x)}
listing the times of all identified peaks in phase with the
first peak.
}
}

A list of the arguments of \code{get_peak_times()} is included as an
attribute.
}
\description{
\code{get_peak_times()} locates peaks in periodic, equally spaced time
series with known period. It applies a central moving average to the
raw time series, identifies when peaks occur in the resulting smooth
time series, and subsets those peaks in phase with the first peak.
}
\section{Bandwidth tuning}{

The bandwidths \code{bw_mavg} and \code{bw_peakid} must be iteratively tuned
so as to disqualify peaks caused by noise in \code{x} (choose \code{bw_mavg}
large enough) without disqualifying true peaks (choose \code{bw_peakid}
not too large). Tuning is not done internally.
}

\section{Missed peaks}{

\code{get_peak_times()} will not detect peaks at the edges of \code{x}.
Specifically, it will not detect a peak at index \code{i} if:
\itemize{
\item \code{i <= floor(bw_mavg) + floor(bw_peakid)}
\item \code{i >= length(x) - floor(bw_mavg) - floor(bw_peakid) + 1}
}

\code{get_peak_times()} will miss peaks in other, typically rare
situations. For example, a peak at index \code{i} will be missed
if, by chance, \code{x_mavg[i] = x_mavg[i+1]} (i.e., there are two
points of equal height at the top of the peak instead of a
single highest point). See the definition of \code{bw_peakid}.
}

\examples{
# Simulate 20 years of disease incidence,
# observed weekly
par_list <- make_par_list(dt_weeks = 1)
df <- make_data(
  par_list = par_list,
  n = 20 * 365 / 7, # 20 years is ~1042 weeks
  with_dem_stoch = TRUE,
  seed = 5
)

# Plot incidence time series, and note the
# apparent 1-year period
plot(Z ~ t_years, df,
  type = "l",
  xlab = "Time (years)",
  ylab = "Incidence"
)

# Find peaks in incidence time series
peaks <- get_peak_times(
  x = df$Z,
  period = 365 / 7, # 1 year is ~52 weeks
  bw_mavg = 6,
  bw_peakid = 8
)

# Verify that peaks were identified correctly
abline(v = df$t_years[peaks$all], lty = 2, col = "red")

}
\references{
deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
time-varying transmission rates for infectious diseases.
}
