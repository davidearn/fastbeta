#' \loadmathjax
#' Deconvolve a time series of reported incidence
#'
#' @description
#' Deconvolves a reported incidence time series using a Richardson-Lucy
#' iteration, generating an estimate of the underlying incidence time series.
#'
#' @param x Reported incidence time series. A numeric vector giving the
#'   number of cases reported on each day from day `0` to day `length(x)-1`.
#' @param delay_dist Distribution of the number of days from infection to
#'   reporting. A numeric vector. `delay_dist[i]` is the probability that
#'   an infection that is eventually reported is reported after `i-1` days.
#'   An error is thrown if `sum(delay_dist) != 1`.
#' @param p Reporting probability. Either
#'   (i) a numeric scalar giving a constant positive probability that an
#'   infection is reported, or
#'   (ii) a numeric vector giving a positive probability for each day
#'   from day `-b` to day `length(x)-1`,
#'   where `b = max(which(delay_dist > 0)) - 1` is the maximum number of
#'   days from infection to reporting. In the latter case, `p[i]` is the
#'   probability that an infection on day `i-1-b` is reported.
#' @param it_max Number of iterations. A non-negative integer.
#'
#' @return
#' A list with elements:
#'
#' \describe{
#'   \item{`times`}{A numeric vector equal to `-b:(length(x)-1)`,
#'     where `b = max(which(delay_dist > 0)) - 1` is the maximum number
#'     of days from infection to reporting, giving time points (in days)
#'     for deconvolved incidence.
#'   }
#'   \item{`poisson`}{A numeric matrix. `poisson[i, j]` is the estimate
#'     after `j-1` iterations of the expected number of infections on day
#'     `times[i]` that are eventually reported. See Details for the initial
#'     estimate.
#'   }
#'   \item{`inc`}{A numeric matrix. `inc[i, j]` is the expected number
#'     of infections on day `times[i]` conditional on `poisson[, j]`.
#'     `inc[, j]` gives deconvolved incidence after `j-1` iterations.
#'   }
#'   \item{`inc_rep`}{A numeric matrix. `inc_rep[i, j]` is the expected
#'     number of cases reported on day `times[i]` conditional on `poisson[, j].`
#'   }
#'   \item{`chi2`}{A numeric vector. `chi2[j]` is the normalized chi-squared
#'     statistic corresponding to `inc_rep[, j]`, computed as
#'     `mean((x - y)^2 / y)`, where `y = inc_rep[b+1:length(x), j]`.
#'   }
#'   \item{`it_chi2_lt1`}{An integer scalar equal to `min(which(chi2 < 1)) - 1`
#'     if `any(chi2 < 1)` is `TRUE` and `NA` otherwise, giving the number of
#'     iterations performed before the chi-squared criterion suggested
#'     by \insertCite{Gold+09;textual}{fastbeta} (namely `chi2 < 1`) was
#'     satisfied. Note that the corresponding incidence time series is
#'     `inc[, 1+it_chi2_lt1]`.
#'   }
#'   \item{`x_pad`}{A numeric vector equal to `c(rep(NA, b), x)`.
#'     A convenience allowing one to plot the observed reported incidence
#'     time series with `lines(times, x_pad, ...)`.
#'   }
#' }
#'
#' Note that `poisson`, `inc`, and `inc_rep` all have `b+length(x)` rows
#' (corresponding to the elements of `times`) and `1+it_max` columns
#' (one for each iteration, starting from an initial value).
#'
#' The list has attributes `call` and `arg_list`, making it reproducible with
#' `eval(call)` or `do.call(deconvolve, arg_list)`.
#'
#' @details
#' # Details
#'
#' ## 1. Differences compared to \insertCite{Gold+09;textual}{fastbeta}
#' The algorithm follows \insertCite{Gold+09;textual}{fastbeta} with the
#' following differences:
#' * Reported incidence is observed from day 0 to day `length(x)-1`,
#'   not from day 1 to day `length(x)`.
#' * Infection and reporting can happen on the same day. A delay
#'   of zero days can be assigned nonzero probability by setting
#'   `delay_dist[1] > 0`.
#' * The initial value for the Richardson-Lucy iteration is obtained by
#'   binning all `sum(x)` cases between day `-b` and day `length(x)-1`,
#'   where `b = max(which(delay_dist > 0)) - 1` is the maximum number
#'   of days from infection to reporting. Binning follows `delay_dist`,
#'   so that `delay_dist[i]` is the fraction of `x[j]` placed in the
#'   bin for day `j-i`. This differs from
#'   \insertCite{Gold+09;textual}{fastbeta}, whose initial value is `x`
#'   shifted backwards in time by `s = which.max(delay_dist) - 1` days.
#'   An undesired consequence of their choice is that deconvolved
#'   incidence is necessarily zero from day `-b` to day `-s-1` and
#'   from day `length(x)-s` to day `length(x)-1`.
#'
#' ## 2. Noise in deconvolved incidence and choosing when to stop
#' \insertCite{Gold+09;textual}{fastbeta} suggest to stop iterating
#' once `chi2 < 1` and to keep the last generated incidence time series
#' (see Value, under `chi2` and `it_chi2_lt1`). The reason for this
#' criterion is that deconvolved incidence tends to acquire undesired
#' noise after several iterations, and it is typically better to stop
#' before this happens;
#' see the supplement to \insertCite{Gold+09;textual}{fastbeta} for
#' a discussion. However, this criterion is not guaranteed to produce
#' an optimal result, hence a decision about which incidence time series
#' to keep should be based on a graphical exploration of the entire
#' `deconvolve()` output (see Examples).
#'
#' ## 3. Missing values in output
#' The last `a` rows of `poisson` and `inc` are `NA`, where
#' `a = min(which(delay_dist > 0)) - 1` is the minimum number of days from
#' infection to reporting. The probability that an infection between day
#' `length(x)-a` and day `length(x)-1` is reported between day `0` and day
#' `length(x)-1` is zero, so the supplied reported incidence time series
#' provides no information about incidence on those `a` days.
#'
#' ## 4. Time step
#' The algorithm does not explicitly require a time step of one day. Longer
#' time steps are tolerated provided `x`, `delay_dist`, and `p` are defined
#' appropriately.
#'
#' @references
#' \insertRef{Gold+09}{fastbeta}
#'
#' @examples
#' x <- 400 * exp(-seq(-2, 2, by = 0.1)^2)
#' delay_dist <- rep(1 / 11, 11)
#' p <- 0.5
#' it_max <- 20
#' dcv_out <- deconvolve(x, delay_dist, p, it_max)
#'
#' plot.new()
#' plot.window(xlim = c(-10, 40), ylim = c(0, 1000))
#' box()
#' axis(side = 1, at = seq(-10, 40, by = 10))
#' title(xlab = "Time (days)")
#' axis(side = 2, at = seq(0, 1000, by = 200), las = 1)
#' title(ylab = "Deconvolved incidence")
#' with(dcv_out, {
#'   for (j in 1:(1+it_max)) {
#'     # Deconvolved incidence after `j` iterations
#'     lines(times, inc[, j])
#'   }
#'   # Deconvolved incidence after `it_max` iterations
#'   lines(times, inc[, 1+it_max], lwd = 2, col = "red")
#'   # Deconvolved incidence after `it_chi2_lt1` iterations
#'   lines(times, inc[, 1+it_chi2_lt1], lwd = 2, col = "blue")
#'   # Observed reported incidence
#'   lines(times, x_pad, lwd = 2, col = "green")
#' })
#'
#' plot.new()
#' plot.window(xlim = c(-10, 40), ylim = c(0, 500))
#' box()
#' axis(side = 1, at = seq(-10, 40, by = 10))
#' title(xlab = "Time (days)")
#' axis(side = 2, at = seq(0, 500, by = 100), las = 1)
#' title(ylab = "Expected reported incidence")
#' with(dcv_out, {
#'   for (j in 1:(1+it_max)) {
#'     # Expected reported incidence after `j` iterations
#'     lines(times, inc_rep[, j])
#'   }
#'   # Expected reported incidence after `it_max` iterations
#'   lines(times, inc_rep[, 1+it_max], lwd = 2, col = "red")
#'   # Expected reported incidence after `it_chi2_lt1` iterations
#'   lines(times, inc_rep[, 1+it_chi2_lt1], lwd = 2, col = "blue")
#'   # Observed reported incidence
#'   lines(times, x_pad, lwd = 2, col = "green")
#' })
#'
#' # Normalized chi-squared statistic
#' plot(0:it_max, dcv_out$chi2, las = 1, xlab = "Iterations", ylab = "chi2")
#' abline(h = 1, lty = 2)
#'
#' @export
deconvolve <- function(x,
                       delay_dist = c(1),
                       p = 1,
                       it_max = 20L) {
  # Stop if `delay_dist` is not valid
  if (sum(delay_dist) != 1) {
    stop("`delay_dist` must sum to 1.")
  }

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Cases are reported at times `0:(N-1)`
  N <- length(x)

  # Time from infection to reporting is at most `b` days
  b <- max(which(delay_dist > 0)) - 1

  # Pad `x` and `delay_dist` with zeros to length `b + N`
  # (convenient for matrix operations)
  x <- c(rep(0, b), x)
  delay_dist <- c(delay_dist[1:(b+1)], rep(0, N - 1))

  # Create a lower triangular transition matrix satisfying
  # `L[i,j] == delay_dist[i-j+1]`
  L <- array(rep(c(delay_dist, 0), b + N), dim = rep(b + N, 2))
  L[upper.tri(L, diag = FALSE)] <- 0

  # Create an upper triangular transition matrix satisfying
  # `U[i,j] == delay_dist[j-i+1]`
  U <- t(L)

  # Compute the probability `q[i]` that an infection on day `i-b-1`
  # is reported on day `j` for some `j` in `0:(N-1)`
  v <- c(rep(0, b), rep(1, N))
  q <- U %*% v
  qnz <- q > 0

  # Preallocate space for all estimates of the Poisson parameter vector
  # and the corresponding time series of (expected) reported incidence
  P <- matrix(NA, nrow = b + N, ncol = 1 + it_max)
  E <- P

  # Initial estimate of the Poisson parameter vector is `x` binned
  # according to `delay_dist`
  P[qnz, 1] <- U[qnz, ] %*% x
  E[, 1] <- L[, qnz] %*% P[qnz, 1]

  # Carry out the Richardson-Lucy iteration
  for (j in 1:it_max) {
    enz <- E[, j] > 0
    P[qnz, j+1] <- (P[qnz, j] / q[qnz]) *
      (U[qnz, enz] %*% (x[enz] / E[enz, j]))
    E[, j+1] <- L[, qnz] %*% P[qnz, j+1]
  }

  # Compute the normalized chi-squared statistic associated with
  # each time series of (expected) reported incidence
  compute_chi2 <- function(y) mean(((x - y)^2 / y)[b+1:N])
  chi2 <- apply(E, 2, compute_chi2)

  # Recover incidence from each Poisson parameter vector
  Z <- P / p

  out <- list(
    times   = -b:(N-1),
    poisson = P,
    inc_rep = E,
    inc     = Z,
    chi2    = chi2,
    x_pad   = c(rep(NA, b), x[b+1:N]),
    it_chi2_lt1 = if (any(chi2 < 1)) min(which(chi2 < 1)) - 1 else NA
  )
  structure(out, call = match.call(), arg_list = arg_list)
}
