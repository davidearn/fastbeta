#' \loadmathjax
#' Bootstrap CIs on transmission rate estimates
#'
#' @description
#' Constructs bootstrap 95% confidence intervals on [fastbeta()] estimates
#' of a time-varying transmission rate \mjseqn{\beta(t)} (or alternatively
#' on smooth loess fits to potentially noisy [fastbeta()] estimates).
#'
#' @details
#' The raw [fastbeta()] estimate of \mjseqn{\beta(t)} (or a smooth loess
#' fit to the raw estimate) is taken to be the true transmission rate for
#' the purpose of simulating new time series data (using [make_data()])
#' and generating bootstrap estimates of \mjseqn{\beta(t)} from those data.
#' The bootstrap 95% CI is then defined as the 2.5-97.5th percentile interval
#' of the bootstrap estimates.
#'
#' To be precise, if a loess object `y` is not specified in the function
#' call, then the data-generating transmission rate is the linear interpolant
#' of the raw [fastbeta()] estimate given in `x$out`, and the bootstrap
#' estimates are the raw [fastbeta()] estimates obtained from the simulated
#' data (one per simulation).
#'
#' On the other hand, if `y` *is* specified, the data-generating transmission
#' rate is the loess fit to the to the raw [fastbeta()] estimate, and the
#' bootstrap estimates are the loess fits to the raw [fastbeta()] estimates
#' obtained from the simulated data, evaluated at the same time points.
#' Values for all loess smoothing parameters are copied from `y`.
#'
#' Arguments `p` and `delay_dist` can be used to define a model for
#' observation error in the data-generating process. Their default
#' values are such that all infections are reported with no delay
#' between infection and reporting, i.e., such that there is no
#' observation error. In each simulation, deconvolution (see [deconvol()])
#' is employed to reconstruct incidence from reported incidence,
#' then the data-generating transmission rate is estimated from
#' deconvolved incidence. The iterative deconvolution algorithm is
#' stopped when the chi-squared criterion is satisfied or when the
#' prescribed maximum number of iterations (25) is reached
#' (whichever occurs first). For details, see Value of [deconvol()].
#'
#' @param x A fastbeta object defining a [fastbeta()] estimate
#'   of \mjseqn{\beta(t)}. Must satisfy `x$method %in% c("si", "sei")`.
#' @param y An optional loess object defining a loess fit to the
#'   [fastbeta()] estimate of \mjseqn{\beta(t)} (see Details).
#'   [try_loess()] can help define a reasonable value for `y`
#'   (see Examples).
#' @param n An integer scalar. The desired number of bootstrap
#'   estimates of \mjseqn{\beta(t)}.
#' @param p A numeric vector of length `nrow(x$out)`. `p[i]` is the
#'   positive probability that an infection between times `x$out$t[i-1]`
#'   and `x$out$t[i]` is eventually reported. Alternatively, a numeric
#'   scalar indicating a constant positive probability.
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
#'     bootstrap estimates of \mjseqn{\beta(t)}. Equal to `x$out$t`.
#'   }
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
#'     Equal to `ncol(mat)`.
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
# Fit a loess curve to the `fastbeta()` estimate.
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
#' # Construct a bootstrap 95% confidence interval
#' # on the loess estimate, conditional on a given
#' # observation model
#' bsbeta_out <- bsbeta(fastbeta_out, my_loess,
#'   n = 100,
#'   p = 0.5, # reporting probability: 50%
#'   delay_dist = c(0, 0.5, 0.5) # reporting delay: 1 dt or 2 dt
#'                               # with equal probability
#' )
#' plot(bsbeta_out)
#'
#' @seealso [methods for class "bsbeta"][bsbeta-methods],
#'   [fastbeta()], [try_loess()], [make_data()]
#'
#' @export
#' @importFrom stats approxfun quantile sd predict
bsbeta <- function(x, y = NULL, n = 100, p = 1, delay_dist = c(1)) {
  if (missing(x)) {
    stop("Missing argument `x`.")
  } else if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  } else if (!x$method %in% c("si", "sei")) {
    stop("`x$method` must be `\"si\"` or `\"sei\"`.")
  }
  if (!is.null(y) && !inherits(y, "loess")) {
    stop("`y` must be a loess object or `NULL`.")
  }
  if (!is.numeric(n) || length(n) != 1 || !isTRUE(n >= 1)) {
    stop("`n` must be a positive integer scalar.")
  }
  if (!is.numeric(p) || !length(p) %in% c(1, nrow(x$out)) ||
        !isTRUE(all(p > 0))) {
    stop("`p` must be a positive numeric vector of length 1 or `nrow(x$out)`.")
  }
  if (!is.numeric(delay_dist) || length(delay_dist) < 1 ||
        !isTRUE(all(delay_dist >= 0))) {
    stop("`delay_dist` must be a non-negative numeric vector.")
  } else if (!any(delay_dist > 0)) {
    stop("`delay_dist` must have at least one positive element.")
  }

  ## Setup ---------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Reserve variable name `n` for `make_data()` argument
  n_bootstrap <- floor(n)


  ## Construct a valid call to `make_data()` using ... -------------------
  ## information stored in the fastbeta object `x`
  ## and loess object `y`

  # List of parameter values
  # * `make_data()` needs additional elements `dt_days` and `N0`.
  # * The precise values have no effect on the `bsbeta()` output,
  #   so dummy values are chosen.
  par_list <- x$arg_list$par_list
  par_list$dt_days <- 1
  par_list$N0 <- with(par_list,
    switch(x$method, si = S0 + I0, sei = S0 + E0 + I0)
  )

  # Time series length
  n <- nrow(x$out)

  # Whether to include demographic stochasticity
  with_ds <- TRUE

  # Whether to use an SIR or SEIR model
  model <- switch(x$method, si = "sir", sei = "seir")

  # Per capita natural mortality rate:
  # * Defined as the linear interpolant of the
  #   PCNMR time series in `x`.
  f <- approxfun(x = x$out$t, y = x$out$mu, rule = 2)
  mu <- function(s, par_list) f(s)

  # Birth rate:
  # * Defined as the linear interpolant of a 2-point
  #   moving average applied the births time series
  #   in `x`.
  # * `B[i] + B[i+1]` is the number of births in the
  #   two observation intervals between times `t[i-1]`
  #   and `t[i+1]`, so `(B[i] + B[i+1]) / 2` gives a
  #   centered estimate of the birth rate per interval
  #   at time `t[i]`.
  gy <- with(x$out, (B + c(B[-1], B[n])) / 2)
  g <- approxfun(x = x$out$t, y = gy, rule = 2)
  nu <- function(s, par_list) g(s)

  # Transmission rate:
  # * Defined as the linear interpolant of the
  #   time series in `x` if `y` was not supplied.
  # * Defined as the loess fit in `y` if `y` was
  #   supplied.
  if (is.null(y)) {
    h <- approxfun(x = x$out$t, y = x$out$beta, rule = 2)
    beta <- function(s, par_list) h(s)
  } else {
    beta <- function(s, par_list) predict(y, s)
  }

  # Reporting probability:
  # * A vector of length `n` is needed even if
  #   the probability is constant.
  if (length(p) == 1) {
    p <- rep(p, n)
  }

  # Reporting delay distribution:
  # * The vector must sum to 1.
  delay_dist <- delay_dist / sum(delay_dist)

  # Maximum reporting delay
  b <- max(which(delay_dist > 0)) - 1


  ## Generate `n_bootstrap` bootstrap estimates ... ----------------------
  ## of the transmission rate

  # Matrix of bootstrap estimates
  mat <- replicate(n_bootstrap, {

    # Data frame containing simulated time series data
    data <- make_data(par_list, n, with_ds, model, mu, nu, beta, p, delay_dist)

    # Data frame containing a deconvolved incidence time series
    deconvol_out <- deconvol(data$C,
      delay_dist = delay_dist,
      p = c(rep(NA, b), p),
      it_max = 25,
      simple = TRUE
    )

    # Data frame to be passed to `fastbeta()`
    data <- data[c("t", "B", "mu")]
    data$Z <- deconvol_out$inc[b+1:nrow(data)]

    # fastbeta object containing raw transmission rate estimate
    fastbeta_out <- fastbeta(data, par_list, method = x$method)

    # Bootstrap estimate of the transmission rate
    if (is.null(y)) {
      fastbeta_out$out$beta
    } else {
      loess_out <- loess(beta ~ t,
        data    = fastbeta_out$out,
        span    = y$pars$span,
        degree  = y$pars$degree,
        control = loess.control(surface = "direct")
      )
      predict(loess_out, x$out$t)
    }
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
