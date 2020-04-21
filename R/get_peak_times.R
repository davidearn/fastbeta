#' Find times of peaks in periodic time series
#'
#' `get_peak_times()` locates peaks in periodic, equally spaced time
#' series with known period. It applies a central moving average to the
#' raw time series, identifies when peaks occur in the resulting smooth
#' time series, and subsets those peaks in phase with the first peak.
#'
#' @section Bandwidth tuning:
#' The bandwidths `bw_mavg` and `bw_peakid` must be iteratively tuned
#' so as to disqualify peaks caused by noise in `x` (choose `bw_mavg`
#' large enough) without disqualifying true peaks (choose `bw_peakid`
#' not too large). Tuning is not done internally.
#'
#' @section Missed peaks:
#' `get_peak_times()` will not detect peaks at the edges of `x`.
#' Specifically, it will not detect a peak at index `i` if:
#'
#' * `i <= floor(bw_mavg) + floor(bw_peakid)`
#' * `i >= length(x) - floor(bw_mavg) - floor(bw_peakid) + 1`
#'
#' `get_peak_times()` will miss peaks in other, typically rare
#' situations. For example, a peak at index `i` will be missed
#' if, by chance, `x_mavg[i] = x_mavg[i+1]` (i.e., there are two
#' points of equal height at the top of the peak instead of a
#' single highest point). See the definition of `bw_peakid`.
#'
#' @param x Numeric vector. A periodic, equally spaced time series.
#'   Time points are taken to be `seq_along(x)`.
#' @param period Numeric scalar. The period of `x` in units of the
#'   observation interval.
#' @param bw_mavg Numeric scalar. Bandwidth for the moving average
#'   applied to `x`. `x_mavg[i]` will be the mean of observations
#'   `x[j]` with
#'   \ifelse{latex}{\out{$|i - j| \leq bandwidth$}}{\ifelse{html}{\out{&vert;<i>i</i> &minus; <i>j</i>&vert; &le; bandwidth}}{|i - j| <= bandwidth}}.
#' @param bw_peakid Numeric scalar. Bandwidth for peak identification.
#'   `i` will be a peak time if and only if `x_mavg[i] > x_mavg[j]`
#'   for all `j` with
#'   \ifelse{latex}{\out{$0 < |i - j| \leq bandwidth$}}{\ifelse{html}{\out{0 &lt; &vert;<i>i</i> &minus; <i>j</i>&vert; &le; bandwidth}}{0 < |i - j| <= bandwidth}}.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{`all`}{Numeric vector. A subset of `seq_along(x)`
#'     listing the times of all identified peaks.
#'   }
#'   \item{`phase`}{Numeric vector. A subset of `seq_along(x)`
#'     listing the times of all identified peaks in phase with the
#'     first peak.
#'   }
#' }
#'
#' A list of the arguments of `get_peak_times()` is included as an
#' attribute.
#'
#' @examples
#' # Simulate 20 years of disease incidence,
#' # observed weekly
#' par_list <- make_par_list(dt_weeks = 1)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7, # 20 years is ~1042 weeks
#'   with_dem_stoch = TRUE,
#'   seed = 5
#' )
#'
#' # Plot incidence time series, and note the
#' # apparent 1-year period
#' plot(Z ~ t_years, df,
#'   type = "l",
#'   xlab = "Time (years)",
#'   ylab = "Incidence"
#' )
#'
#' # Find peaks in incidence time series
#' peaks <- get_peak_times(
#'   x = df$Z,
#'   period = 365 / 7, # 1 year is ~52 weeks
#'   bw_mavg = 6,
#'   bw_peakid = 8
#' )
#'
#' # Verify that peaks were identified correctly
#' abline(v = df$t_years[peaks$all], lty = 2, col = "red")
#'
#' @references
#' deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
#' time-varying transmission rates for infectious diseases.
#' 
#' @md
#' @export
get_peak_times <- function(x, period, bw_mavg, bw_peakid) {

## 1. Apply moving average ---------------------------------------------

# Fractional part of bandwidth can be discarded
# without affecting number of observations in band
bw_mavg <- floor(bw_mavg)

# Number of observations in moving average band
m <- 2 * bw_mavg + 1

# `filter()` applies a moving average to `x`. It
# returns a vector of equal length containing the
# resulting smooth time series.
x_mavg <- as.vector(
  stats::filter(x,
    filter = rep(1 / m, m),
    method = "convolution",
    sides  = 2
  )
)


## 2. Locate all peaks -------------------------------------------------

# Fractional part of bandwidth can be discarded
# without affecting number of observations in band
bw_peakid <- floor(bw_peakid)

# Number of observations in peak identification band
m <- 2 * bw_peakid + 1

# Pads ensure that `embed()` returns incomplete
# bands, e.g., the band centered at `x_mavg[1]`
pad <- rep(NA, bw_peakid)
x_mavg_padded <- c(pad, x_mavg, pad)

# `embed()` applied to `x_mavg_padded` returns a
# matrix whose `i`th row is the band centered at
# `x_mavg[i]`
bands <- stats::embed(x_mavg_padded, m)

# Central element of each band
bands_mid <- bands[, bw_peakid+1]
# Remaining elements
bands_sides <- bands[, -(bw_peakid+1)]
# Maximum of remaining elements
bands_sides_max <- apply(bands_sides, 1, max)

# Indices of bands in which the central element exceeds
# the remaining elements. These are defined as peak times.
peak_times_all <- which(bands_mid > bands_sides_max)


## 3. Subset peaks in phase with first peak ----------------------------

# First peak time
tA <- peak_times_all[1]
# Last observation time
tf <- length(x)
# Number of cycles in between
num_cycles <- floor((tf - tA) / period)

# Peak times in phase with first peak time, in the
# hypothetical case of exactly periodic dynamics
peak_times_phase <- tA + (0:num_cycles) * period

# Dynamics in the time series `x_mavg` may not be
# exactly periodic. Replace each hypothetical time
# in `peak_times_phase` with the nearest observed
# time in `peak_times_all`.
peak_times_phase <- sapply(peak_times_phase,
  function(x) peak_times_all[which.min(abs(peak_times_all - x))]
)


out <- list(
  all   = peak_times_all,
  phase = peak_times_phase
)
attr(out, "arg_list") <-
  as.list(environment())[names(formals(get_peak_times))]
out
}
