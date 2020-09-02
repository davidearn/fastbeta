#' \loadmathjax
#' Bootstrap CIs on transmission rate estimates
#'
#' @description
#' Constructs bootstrap 95% confidence intervals on [fastbeta()] estimates
#' of a time-varying transmission rate \mjseqn{\beta(t)} (or alternatively
#' on loess fits to [fastbeta()] estimates).
#'
#' @details
#' The [fastbeta()] estimate of \mjseqn{\beta(t)} (or a loess fit to the
#' estimate) is taken to be the true transmission rate for the purpose of
#' simulating new time series data and generating bootstrap estimates of
#' \mjseqn{\beta(t)} from those data. The bootstrap 95% CI is then defined
#' as the 2.5-97.5th percentile interval of the bootstrap estimates.
#'
#' Arguments `p` and `delay_dist` can be used to define a model for
#' observation error in the data-generating process. Their default
#' values are such that all infections are reported with no delay
#' between infection and reporting, i.e., such that there is no
#' observation error.
#'
#' @param x A fastbeta object defining a [fastbeta()] estimate
#'   of \mjseqn{\beta(t)}.
#' @param y An optional loess object defining a loess fit
#'   to the [fastbeta()] estimate (see Details). [try_loess()]
#'   can help define a reasonable value for `y` (see Examples).
#' @param n An integer scalar. The desired number of bootstrap
#'   estimates of \mjseqn{\beta(t)}.
#' @param p A numeric vector of length `nrow(x$out)-1`. `p[i]` is the
#'   probability that an infection between times `x$out$t[i]` and
#'   `x$out$t[i+1]` is eventually reported. Alternatively, a numeric
#'   scalar indicating a constant probability.
#' @param delay_dist A numeric vector. The distribution of the integer
#'   number of observation intervals between infection and reporting.
#'   `delay_dist[i]` is the probability that an infection that is
#'   eventually reported is reported after `i-1` observation intervals.
#'   `delay_dist` is replaced with `delay_dist / sum(delay_dist)` in
#'   the event that `sum(delay_dist) != 1`.
#'
#' @return
#' A bsbeta object. A list with elements:
#'
#' \describe{
#'   \item{`times`}{A numeric vector listing time points for the
#'     bootstrap estimates of \mjseqn{\beta(t)}}. Equal to `x$out$t`.
#'   \item{`beta`}{A numeric vector. If `y = NULL`, the [fastbeta()]
#'     estimate of \mjseqn{\beta(t)}, equal to `x$out$beta`. Otherwise,
#'     the loess fit evaluated at `times`, equal to `predict(y, times)`.
#'   }
#'   \item{`mat`}{A numeric matrix whose `j`th column is the
#'     `j`th bootstrap estimate of \mjseqn{\beta(t)}.
#'   }
#'   \item{`ci95`}{A numeric matrix with 2 columns giving the lower
#'     and upper endpoints of the bootstrap 95% confidence interval.
#'     Equal to
#'     `t(apply(mat, 1, function(x) if (all(is.finite(x))) quantile(x, probs = c(0.025, 0.975)) else rep(NA, 2)))`.
#'   }
#'   \item{`mean`}{A numeric vector. The mean bootstrap estimate.
#'     Equal to `rowMeans(mat)`.
#'   }
#'   \item{`sd`}{A numeric vector. The standard deviation of
#'     the bootstrap estimates. Equal to `apply(mat, 1, sd)`.
#'   }
#'   \item{`n`}{An integer scalar. The number of bootstrap estimates.
#'     Equal to `n` in the function call and also `ncol(mat)`.
#'   }
#'   \item{`call`}{The function call. The bsbeta object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call.
#'     The bsbeta object is reproducible with `do.call(bsbeta, arg_list)`.
#'   }
#' }
#'
#' @examples
#' # Simulate time series data using an SIR model
#' # with seasonally forced transmission rate
#' pl <- make_par_list(model = "sir")
#' df <- make_data(pl, with_ds = TRUE, model = "sir")
#'
#' # Estimate the seasonally forced transmission rate
#' # using the SI method
#' fastbeta_out <- fastbeta(df, pl, method = "si")
#' plot(fastbeta_out)
#'
#' # Fit a loess curve to the `fastbeta()` estimate.
#' # Try different values for the smoothing parameter
#' # and choose one that produces a reasonable fit.
#' try_loess_out <- try_loess(beta ~ t,
#'   data = fastbeta_out$out,
#'   q = seq(20, 90, by = 10),
#'   control = loess.control(surface = "direct"),
#'   plot = TRUE
#' )
#' my_loess <- try_loess_out[["q50"]] # `q = 50` seems reasonable
#'
#' # Construct bootstrap 95% confidence intervals
#' # on the loess estimate
#' bsbeta_out <- bsbeta(fastbeta_out, my_loess, n = 100)
#' plot(bsbeta_out)
#'
#' @export
#' @importFrom stats approxfun quantile sd predict
bsbeta <- function(x, y = NULL, n = 100, p = 1, delay_dist = c(1)) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  if (!is.null(y) && !inherits(y, "loess")) {
    stop("`y` must be a loess object or `NULL`.")
  }
  arg_list <- as.list(environment())
  n_bootstrap <- n
  par_list <- x$arg_list$par_list
  par_list$dt_days <- NA
  par_list$N0 <- with(par_list,
    switch(x$method, si = S0 + I0, sei = S0 + E0 + I0)
  )
  n <- nrow(x$out)
  with_ds <- TRUE
  model <- switch(x$method, si = "sir", sei = "seir")
  f <- approxfun(x = x$out$t, y = x$out$mu, rule = 2)
  mu <- function(s, par_list) f(s)
  g <- approxfun(x = x$out$t, y = with(x$out, (B + c(B[-1], B[n])) / 2), rule = 2)
  nu <- function(s, par_list) g(s)
  if (is.null(y)) {
    h <- approxfun(x = x$out$t, y = x$out$beta, rule = 2)
    beta <- function(s, par_list) h(s)
  } else {
    beta <- function(s, par_list) predict(y, s)
  }
  if (length(p) != n - 1) {
    pconst <- p[1]
    p <- rep(pconst, n - 1)
  }
  delay_dist <- delay_dist / sum(delay_dist)
  mat <- replicate(n_bootstrap, {
    data <- make_data(par_list, n, with_ds, model, mu, nu, beta, p, delay_dist)
    fastbeta(data, par_list, method = x$method)$out$beta
  })
  out <- list(
    times = x$out$t,
    beta  = if (is.null(y)) x$out$beta else predict(y, x$out$t),
    mat   = mat,
    ci95  = t(apply(mat, 1, function(x) if (all(is.finite(x))) quantile(x, probs = c(0.025, 0.975)) else rep(NA, 2))),
    mean  = rowMeans(mat),
    sd    = apply(mat, 1, sd),
    n     = ncol(mat),
    call  = match.call(),
    arg_list = arg_list
  )
  colnames(out$ci95) <- c("lower", "upper")
  structure(out, class = c("bsbeta", "list"))
}
