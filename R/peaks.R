#' \loadmathjax
#' Locate peaks in time series
#'
#' @description
#' Locates peaks in time series with observations at equally spaced
#' time points \mjseqn{t_i = t_0 + i \Delta t}. For roughly periodic
#' time series, identifies those peaks in phase with the first peak.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#' The aim is to index peaks in the supplied time series
#' \mjseqn{\lbrace (t_i,x_i) \rbrace_{i=0}^{n-1}}. To ensure that
#' this output is robust to noise, a \mjseqn{(2 \ell_1 + 1)}-point
#' central moving average is applied, generating a smoother time series
#' \mjseqn{\out{\lbrace (t_i,\bar{x}_i) \rbrace_{i=\ell_1}^{n-1-\ell_1}}}
#' with
#'
#' \mjsdeqn{\out{\bar{x}_i = \frac{1}{2 \ell_1 + 1} \sum_{j=-\ell_1}^{\ell_1} x_{i+j}\,.}}
#'
#' Peaks in
#' \mjseqn{\out{\lbrace (t_i,\bar{x}_i) \rbrace_{i=\ell_1}^{n-1-\ell_1}}}
#' are found by taking temporal neighbourhoods containing
#' \mjseqn{(2 \ell_2 + 1)} points and checking whether the central point
#' is a local maximum. Specifically, \mjseqn{(t_k,\bar{x}_k)} is defined
#' as a peak if and only if \mjseqn{k \in \mathcal{I}}, where
#'
#' \mjsdeqn{\out{\mathcal{I} = \Big\lbrace i \in \lbrace\ell_1+\ell_2,\ldots,n-1-\ell_1-\ell_2\rbrace : \bar{x}_i > \bar{x}_{i \pm j}\,\text{for all}\,j = 1,\ldots,\ell_2 \Big\rbrace\,.}}
#'
#' If the unobserved function \mjseqn{x(t)} is \mjseqn{T}-periodic,
#' then a secondary aim is to identify all peaks in phase with the
#' first peak at \mjseqn{\out{(t_a,\bar{x}_a)}},
#' where \mjseqn{a = \min(\mathcal{I})}.
#' \mjseqn{\out{(t_k,\bar{x}_k)}} is defined as a peak in phase with
#' the first peak if and only if
#' \mjseqn{k \in \mathcal{I}_\text{phase} \subset \mathcal{I}}, where
#'
#' \mjsdeqn{\out{\mathcal{I}_\text{phase} = \Big\lbrace \arg\min_{i \in \mathcal{I}} \big|t_i - (t_a + jT)\big| : j = 0,\ldots,m \Big\rbrace\,.}}
#'
#' Here, \mjseqn{m} is the number of complete cycles between times
#' \mjseqn{t_a} and \mjseqn{t_{n-1-\ell_1}}, given by
#'
#' \mjsdeqn{m = \left\lfloor \frac{t_{n-1-\ell_1} - t_a}{T} \right\rfloor\,.}
#'
#' Hence a peak is considered in phase with the first peak if and
#' only if it occurs nearest to an integer multiple of the period
#' after the first peak.
#'
#' ## 2. Missed peaks
#' The procedure outlined in Algorithm will not detect peaks at the edges
#' of \mjseqn{\lbrace (t_i,x_i) \rbrace_{i=0}^{n-1}}. Specifically, it will
#' not detect a peak at time \mjseqn{t_k} if \mjseqn{k < \ell_1+\ell_2} or
#' \mjseqn{k > n-1-\ell_1-\ell_2}. It can miss peaks in other (unlucky)
#' situations. For example, a peak at time \mjseqn{t_k} will be missed if,
#' by chance, \mjseqn{\out{\bar{x}_k = \bar{x}_{k+1}}} (i.e., if, after
#' smoothing, there are two points of equal height at the top of a peak
#' instead of a single highest point).
#'
#' ## 3. Bandwidth tuning
#' The bandwidths `bw1` and `bw1` (denoted by \mjseqn{\ell_1} and
#' \mjseqn{\ell_2} in Algorithm) must be tuned in order to disqualify
#' spurious peaks in \mjseqn{\lbrace (t_i,x_i) \rbrace_{i=0}^{n-1}}
#' caused by noise, without disqualifying true peaks. Disqualifying
#' spurious peaks is a matter of smoothing to a sufficient degree
#' (choosing `bw1` large enough) and including a nontrivial number
#' of points in the temporal neighbourhoods that distinguish local
#' maxima from other points (choosing `bw2` large enough). Not
#' disqualifying true peaks is a matter of not including more than
#' one peak in those temporal neighbourhoods (choosing `bw2` less
#' than `0.5 * num_dt_in_period / num_peaks_in_period`, if peaks are
#' evenly spaced within a period).
#'
#' Reasonable results are typically obtained by:
#' * Choosing the smallest `bw1` such that
#'   \mjseqn{\out{\lbrace (t_k,\bar{x}_k) \rbrace_{k=\ell_1}^{n-1-\ell_1}}}
#'   is smooth around the peaks, which can be found by experimenting with
#'   [stats::filter()] (see Examples).
#' * Choosing `bw2` around
#'   `0.4 * num_dt_in_period / num_peaks_in_period`.
#'
#' @param x \mjseqn{\lbrace\,x_i\,\rbrace}
#'   A numeric vector defining an equally spaced time series.
#' @param bw1 \mjseqn{\lbrace\,\ell_1\,\rbrace}
#'   An integer scalar. Bandwidth for the central moving average
#'   applied to `x`, so that `xbar[i] = mean(x[(i-bw1):(i+bw1)])`
#'   for all `i`.
#' @param bw2 \mjseqn{\lbrace\,\ell_2\,\rbrace}
#'   An integer scalar. Bandwidth for peak definition, so that
#'   a peak occurs at index `i` if and only if
#'   `all(xbar[i] > xbar[c((i-bw2):(i-1), (i+1):(i+bw2))])`.
#' @param period \mjseqn{\lbrace\,T / \Delta t\,\rbrace}
#'   A numeric scalar. Period of `x` in units \mjseqn{\Delta t}.
#'   Necessary only if `x` is roughly periodic and peaks in phase
#'   with the first peak are desired.
#'
#' @return
#' A peaks object. A list with elements:
#'
#' \describe{
#'   \item{`x`}{A numeric vector. The value of `x` in the function call.}
#'   \item{`xbar`}{A numeric vector. The result of applying
#'     a `(2*bw1+1)`-point central moving average to `x`.
#'   }
#'   \item{`all`}{An integer vector. A subset of `seq_along(x)`
#'     listing the index of each peak in `x`.
#'   }
#'   \item{`phase`}{An integer vector. A subset of `all` listing
#'     the index of each peak in phase with the first peak.
#'     `NULL` if `period` is `NULL` in the function call.
#'   }
#'   \item{`call`}{The function call. The peaks object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call.
#'     The peaks object is reproducible with `do.call(peaks, arg_list)`
#'   }
#' }
#'
#' @examples
#' ## Create a noisy 2-periodic time series
#' ## with 2 peaks per period
#' times <- seq(0, 20, by = 0.01)
#' x <- sin(2 * pi * times) + sin(pi * times) + rnorm(times, 0, 0.5)
#' plot(x, type = "l", lwd = 2, col = "grey80", las = 1)
#'
#' ## Smoothing with `bw1 = 20` removes noise
#' ## around the peaks
#' bw1 <- 20
#' m <- 2 * bw1 + 1
#' xbar <- as.vector(
#'   stats::filter(x,
#'     filter = rep(1 / m, m),
#'     method = "convolution",
#'     sides  = 2
#'   )
#' )
#' lines(xbar, lwd = 2)
#'
#' ## Setting `bw2 = 0.4 * num_dt_in_period / num_peaks_in_period`
#' ## tends to yield reasonable results
#' bw2 <- 0.4 * 200 / 2
#'
#' ## Locate peaks in time series
#' peaks_out <- peaks(
#'   x = x,
#'   bw1 = bw1,
#'   bw2 = bw2,
#'   period = 200 # period is 200 observation intervals
#' )
#'
#' ## Verify that peaks were identified correctly:
#' ## * Blue lines are peaks.
#' ## * Pink circles are peaks in phase with first peak.
#' plot(peaks_out)
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' @seealso [methods for class "peaks"][peaks-methods], [ptpi()]
#' @export
#' @importFrom stats filter embed
peaks <- function(x, bw1, bw2, period = NULL) {
  if (missing(x)) {
    stop("Missing argument `x`.")
  } else if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.")
  }
  if (missing(bw1)) {
    stop("Missing argument `bw1`.")
  } else if (missing(bw2)) {
    stop("Missing argument `bw2`.")
  } else if (!is.numeric(bw1) || length(bw1) != 1 || !isTRUE(bw1 >= 0) ||
               !is.numeric(bw2) || length(bw2) != 1 || !isTRUE(bw2 >= 0)) {
    stop("`bw1` and `bw2` must be non-negative numeric scalars.")
  }
  if (!is.null(period) && (!is.numeric(period) || !isTRUE(period > 0))) {
    stop("`period` must be `NULL` or a positive numeric scalar.")
  }

  ### Setup ------------------------------------------------------------

  ## Save arguments in a list
  arg_list <- as.list(environment())

  ## Fractional part of bandwidths can be discarded
  ## without affecting number of observations in band
  bw1 <- floor(bw1)
  bw2 <- floor(bw2)


  ### Apply moving average ---------------------------------------------

  ## Number of observations in moving average band
  m <- 2 * bw1 + 1

  ## `filter()` applies a moving average to `x`. It
  ## returns a vector of equal length containing the
  ## resulting smooth time series.
  xbar <- as.vector(
    filter(x,
      filter = rep(1 / m, m),
      method = "convolution",
      sides  = 2
    )
  )


  ### Locate all peaks -------------------------------------------------

  ## Number of observations in peak definition band
  m <- 2 * bw2 + 1

  ## Pads ensure that `embed()` returns incomplete
  ## bands, e.g., the band centered at `xbar[1]`
  pad <- rep(NA, bw2)
  xbar_padded <- c(pad, xbar, pad)

  ## `embed()` applied to `xbar_padded` returns a
  ## matrix whose `i`th row is the band centered at
  ## `xbar[i]`
  bands <- embed(xbar_padded, m)

  ## Central element of each band
  bands_mid <- bands[, bw2+1]
  ## Remaining elements
  bands_sides <- bands[, -(bw2+1)]
  ## Maximum of remaining elements
  bands_sides_max <- apply(bands_sides, 1, max)

  ## Indices of bands in which the central element exceeds the
  ## remaining elements. These can be considered peak "times".
  peaks_all <- which(bands_mid > bands_sides_max)


  ### Subset peaks in phase with first peak ----------------------------

  if (!is.null(period)) {
    ## First peak time
    ta <- peaks_all[1]
    ## Last observation time
    tf <- length(x) - bw1
    # Number of cycles in between
    num_cycles <- floor((tf - ta) / period)

    ## Peak times in phase with first peak time in the
    ## hypothetical case of exactly periodic dynamics
    peaks_phase <- ta + (0:num_cycles) * period

    ## Dynamics in the time series `xbar` may not be
    ## exactly periodic. Replace each hypothetical time
    ## in `peaks_phase` with the nearest observed time
    ## in `peaks_all`.
    peaks_phase <- sapply(peaks_phase,
      function(x) peaks_all[which.min(abs(peaks_all - x))]
    )
    peaks_phase <- peaks_phase[!duplicated(peaks_phase)]
  }

  out <- list(
    x        = x,
    xbar     = xbar,
    all      = peaks_all,
    phase    = if (is.null(period)) NULL else peaks_phase,
    call     = match.call(),
    arg_list = arg_list
  )
  structure(out, class = c("peaks", "list"))
}
