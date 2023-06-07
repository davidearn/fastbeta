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
  x_mavg = x_mavg,
  all    = peak_times_all,
  phase  = peak_times_phase
)
attr(out, "arg_list") <-
  as.list(environment())[names(formals(get_peak_times))]
out
}
