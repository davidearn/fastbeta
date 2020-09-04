#' Deconvolve a reported incidence time series
#'
#' @description
#' Deconvolves an equally spaced reported incidence time series using a
#' Richardson-Lucy iteration, generating a reconstruction of the underlying
#' incidence time series.
#'
#' @details
#' # Details
#'
#' ## 1. Departures from \insertCite{Gold+09;textual}{fastbeta}
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
#' once `chi2 < 1` and to keep the last generated incidence time
#' series (see `chi2` in Value). The reason for this criterion is that
#' deconvolved incidence can acquire undesired noise after several
#' iterations, and it is typically better to stop before this happens
#' (see the supplement to \insertCite{Gold+09;textual}{fastbeta} for
#' details). However, this criterion is not guaranteed to produce an
#' optimal result, hence a decision about which incidence time series
#' to keep should be based on a graphical exploration of the entire
#' `deconvol()` output (see Examples).
#'
#' ## 3. Missing values in output
#' The last `a` rows of `poisson` and `inc` are `NA`, where
#' `a = min(which(delay_dist > 0)) - 1` is the minimum number of days from
#' infection to reporting. For `a > 0`, the probability that an infection
#' between day `length(x)-a` and day `length(x)-1` is reported between day
#' `0` and day `length(x)-1` is zero, so the supplied reported incidence
#' time series provides no information about incidence on those `a` days.
#'
#' ## 4. Observation interval
#' The algorithm does not explicitly require daily observations.
#' For example, `x` can be described most generally as listing the
#' number of infections reported during each observation interval
#' from interval 0 to interval `length(x)-1`. A different observation
#' interval (e.g., one week instead of one day) can be chosen as long
#' as the choice is consistent across `x`, `delay_dist`, and `p`.
#' The output must be interpreted accordingly.
#'
#' @param x A numeric vector defining an equally spaced reported incidence
#'   time series. Lists the number of cases reported on each day from day 0
#'   to day `length(x)-1`.
#' @param delay_dist A numeric vector. The distribution of the number of days
#'   from infection to reporting. `delay_dist[i]` is the probability that an
#'   infection that is eventually reported is reported after `i-1` days.
#'   `delay_dist` is replaced with `delay_dist / sum(delay_dist)` in the
#'   event that `sum(delay_dist) != 1`.
#' @param p Either
#'   (i) a numeric scalar giving a constant probability that an infection
#'   is reported, or
#'   (ii) a numeric vector giving a probability for each day from day `-b`
#'   to day `length(x)-1`, where `b = max(which(delay_dist > 0)) - 1` is
#'   the maximum number of days from infection to reporting. In the latter
#'   case, `p[i]` is the probability that an infection on day `i-1-b` is
#'   reported. `NA` are tolerated but cause `NA` to appear in the
#'   corresponding rows of `inc` (see `inc` in Value).
#' @param it_max A non-negative integer. The desired number of iterations.
#' @param simple A logical scalar. If `TRUE`, then a data frame with
#'   2 columns is returned instead of a deconvol object. See Value.
#'
#' @return
#' If `simple = FALSE` (the default), a deconvol object.
#' A list with elements:
#'
#' \describe{
#'   \item{`times`}{A numeric vector equal to `-b:(length(x)-1)`,
#'     where `b = max(which(delay_dist > 0)) - 1` is the maximum number
#'     of days from infection to reporting, giving time points in days
#'     for deconvolved incidence.
#'   }
#'   \item{`x_pad`}{A numeric vector equal to `c(rep(NA, b), x)`.
#'     A convenience allowing one to plot the observed reported incidence
#'     time series with `lines(times, x_pad, ...)`.
#'   }
#'   \item{`poisson`}{A numeric matrix. `poisson[i, j]` is the estimate
#'     after `j-1` iterations of the expected number of infections on day
#'     `times[i]` that are eventually reported. See Details 1 for the
#'     initial estimate.
#'   }
#'   \item{`inc`}{A numeric matrix. `inc[i, j]` is the expected number
#'     of infections on day `times[i]` conditional on `poisson[, j]`.
#'     `inc[, j]` gives deconvolved incidence after `j-1` iterations.
#'   }
#'   \item{`inc_rep`}{A numeric matrix. `inc_rep[i, j]` is the expected
#'     number of cases reported on day `times[i]` conditional on `poisson[, j].`
#'   }
#'   \item{`chi2`}{A numeric vector. `chi2[j]` is the normalized
#'     chi-squared statistic corresponding to `inc_rep[, j]`,
#'     computed as
#'     `mean((x - y)^2 / y)`, where `y = inc_rep[b+1:length(x), j]`.
#'     \insertCite{Gold+09;textual}{fastbeta} suggest taking column
#'     `j = min(which(chi2 < 1))` of `inc` to estimate incidence
#'     instead of the last column (see Details 2).
#'   }
#'   \item{`call`}{The function call. The deconvol object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call. The
#'     deconvol object is reproducible with `do.call(deconvol, arg_list)`.
#'   }
#' }
#'
#' `poisson`, `inc`, and `inc_rep` all have `b+length(x)` rows
#' (corresponding to the elements of `times`) and `1+it_max` columns
#' (one for each iteration, starting from an initial value).
#'
#' If `simple = TRUE`, then a data frame whose first column is `times`
#' and whose second column is `inc[, min(which(chi2 < 1))]` or
#' otherwise `inc[, ncol(inc)]`, depending on whether the chi-squared
#' criterion was reached.
#'
#' @references
#' \insertRef{Gold+09}{fastbeta}
#'
#' @examples
#' x <- 400 * exp(-seq(-2, 2, by = 0.1)^2)
#' delay_dist <- rep(1 / 11, 11)
#' p <- 0.5
#' it_max <- 25
#' deconvol_out <- deconvol(x, delay_dist, p, it_max)
#' plot(deconvol_out) # creates 3 plots
#'
#' @seealso [methods for class "deconvol"][deconvol-methods], [convol()]
#' @export
deconvol <- function(x,
                     delay_dist = c(1),
                     p = 1,
                     it_max = 20L,
                     simple = FALSE) {
  ## Save arguments in a list
  arg_list <- as.list(environment())

  ## Normalize `delay_dist`
  delay_dist <- delay_dist / sum(delay_dist)

  ## Cases are reported at times `0:(N-1)`
  N <- length(x)

  ## Time from infection to reporting is at most `b` days
  b <- max(which(delay_dist > 0)) - 1

  ## Pad `x` and `delay_dist` with zeros to length `b + N`
  ## (convenient for matrix operations)
  x <- c(rep(0, b), x)
  delay_dist <- c(delay_dist[1:(b+1)], rep(0, N - 1))

  ## Create a lower triangular transition matrix satisfying
  ## `L[i,j] == delay_dist[i-j+1]`
  L <- array(rep(c(delay_dist, 0), b + N), dim = rep(b + N, 2))
  L[upper.tri(L, diag = FALSE)] <- 0

  ## Create an upper triangular transition matrix satisfying
  ## `U[i,j] == delay_dist[j-i+1]`
  U <- t(L)

  ## Compute the probability `q[i]` that an infection on day `i-b-1`
  ## is reported on day `j` for some `j` in `0:(N-1)`
  v <- c(rep(0, b), rep(1, N))
  q <- U %*% v
  qnz <- q > 0

  ## Preallocate space for all estimates of the Poisson parameter vector
  ## and the corresponding time series of (expected) reported incidence
  P <- matrix(NA, nrow = b + N, ncol = 1 + it_max)
  E <- P

  ## Initial estimate of the Poisson parameter vector is `x`
  ## binned backwards in time according to `delay_dist`
  P[qnz, 1] <- U[qnz, ] %*% x
  E[, 1] <- L[, qnz] %*% P[qnz, 1]

  ## Carry out the Richardson-Lucy iteration
  for (j in 1:it_max) {
    enz <- E[, j] > 0
    P[qnz, j+1] <- (P[qnz, j] / q[qnz]) *
      (U[qnz, enz] %*% (x[enz] / E[enz, j]))
    E[, j+1] <- L[, qnz] %*% P[qnz, j+1]
  }

  ## Compute the normalized chi-squared statistic associated with
  ## each time series of (expected) reported incidence
  compute_chi2 <- function(y) mean(((x - y)^2 / y)[b+1:N])
  chi2 <- apply(E, 2, compute_chi2)

  ## Recover incidence from each Poisson parameter vector
  Z <- P / p

  if (simple) {
    index <- if (any(isTRUE(chi2 < 1))) min(which(chi2 < 1)) else ncol(Z)
    out <- data.frame(times = -b:(N-1), inc = Z[, index])
  } else {
    out <- list(
      times    = -b:(N-1),
      poisson  = P,
      inc_rep  = E,
      inc      = Z,
      chi2     = chi2,
      x_pad    = c(rep(NA, b), x[b+1:N]),
      call     = match.call(),
      arg_list = arg_list
    )
    out <- structure(out, class = c("deconvol", "list"))
  }
  out
}
