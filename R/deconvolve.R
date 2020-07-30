#' \loadmathjax
#' Reconstruct incidence from reported incidence
#'
#' @description
#' `deconvolve()` employs a Richardson-Lucy iteration to deconvolve a reported
#' incidence time series. The supplied delay distribution is used to carry out
#' the deconvolution, and the result is scaled by the supplied
#' reporting probability to generate an incidence time series.
#'
#' @param x Reported incidence time series. A numeric vector giving the
#'   number of cases reported on each day from day `0` to day `length(x)-1`.
#' @param delay_dist Distribution of the number of days from infection to
#'   reporting. A numeric vector. `delay_dist[i]` is the probability that
#'   an infection that is eventually reported is reported after `i-1` days.
#'   An error is thrown if `sum(delay_dist) != 1`.
#' @param p Reporting probability. Either
#'   (i) a numeric scalar giving a positive probability that an infection
#'   is reported, or
#'   (ii) a numeric vector giving a positive probability for each day
#'   from day `-b` to day `length(x)-1`,
#'   where `b = max(which(delay_dist > 0)) - 1` is the maximum number of
#'   days from infection to reporting. In the latter case, `p[i]` is the
#'   probability that an infection on day `i-1-b` is reported.
#' @param max_it Maximum number of iterations. A non-negative integer.
#'
#' @return
#' A numeric vector of length `b+length(x)` giving an estimate of the
#' number of infections on each day from day `-b` to day `length(x)-1`,
#' where `b = max(which(delay_dist > 0)) - 1` is the maximum number of
#' days from infection to reporting. It is the last Richardson-Lucy
#' iterate, scaled by `1 / p`.
#'
#' The returned vector has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(deconvolve, arg_list)`.
#' It has a third attribute `mat`, a matrix containing the full output
#' of the Richardson-Lucy iteration, prior to scaling by `1 / p`. There
#' is one row for each day from day `-b` to day `length(x)-1` and one
#' column for each iteration (not counting the first column, which gives
#' the initial value for the iteration). Hence the returned vector is
#' the last column of `mat`, scaled by `1 / p`.
#'
#' @references
#' \insertRef{Gold+09}{fastbeta}
#'
#' @export
deconvolve <- function(x,
                       delay_dist = c(1),
                       p = rep(1, length(x) + length(delay_dist) - 1),
                       max_it = 20L) {
  # Save arguments in a list
  arg_list <- as.list(environment())
  
  # Cases are reported at times `0:(N-1)`
  N <- length(x)

  # Time from infection to reporting is between `a` and `b`
  delay_dist_support <- which(delay_dist > 0) - 1
  a <- min(delay_dist_support)
  b <- max(delay_dist_support)

  # Time runs from earliest possible infection to latest report
  t_out <- -b:(N-1)

  # Pad `x` and `delay_dist` with zeros to `length(t_out)`. Having
  # parallel vectors will be convenient for matrix operations.
  x <- c(rep(0, b), x)
  delay_dist <- c(delay_dist[1:(b+1)], rep(0, N - 1))

  # Create a lower triangular transition matrix satisfying
  # `L[i,j] == delay_dist[i-j+1]`
  L <- array(rep(c(delay_dist, 0), b + N), dim = rep(b + N, 2))
  L[upper.tri(L, diag = FALSE)] <- 0

  # Create an upper triangular transition matrix satisfying
  # `U[i,j] == delay_dist[j-i+1]`
  U <- t(L)

  # Compute the probability `q[i]` that an infection
  # at time `t_out[i]` is reported at a time in `0:(N-1)`
  v <- c(rep(0, b), rep(1, N))
  q <- U %*% v

  # Retain all estimates of the Poisson parameter vector
  P <- matrix(NA, nrow = b + N, ncol = max_it + 1)

  # Initial estimate is `x` shifted backwards by the expected delay
  delay_mean <- (0:b) * delay_dist[1:(b+1)]
  P[, 1] <- utils::tail(c(x, rep(0, delay_mean)), b + N)

  # Perform Richardson-Lucy iteration until the chi-squared criterion
  # is satisfied or until `max_it` is reached
  chi2 <- 1
  it <- 1
  while (chi2 >= 1 && it <= max_it) {
    x_expected <- L %*% P[, it]
    P[, it+1] <- (P[, it] / q) * (U %*% (x / x_expected))
    chi2 <- sum(utils::tail((x_expected - x)^2 / x_expected, N)) / N
    it <- it + 1
  }

  # Recover an incidence time series from each Poisson parameter vector
  out <- P / p

  attr(out, "call") <- match.call()
  attr(out, "arg_list") <- arg_list
  attr(out, "mat") <- P
  out
}
