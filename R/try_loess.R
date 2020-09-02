#' Explore loess fits to equally spaced time series
#'
#' @description
#' A wrapper of [stats::loess()] for equally spaced time series.
#' Fits a loess curve to the specified time series using each
#' supplied value of the smoothing parameter and optionally plots
#' the fitted curves.
#'
#' @details
#' When applying [stats::loess()] to smooth an equally spaced
#' time series, all possible smoothing kernels are captured by
#' `span = q / nrow(data)` with integers `q`. All other values
#' of `span` (those produced by non-integer `q`) are redundant.
#' See \insertCite{ClevGros+91;textual}{fastbeta} for details.
#'
#' When fitting local quadratic polynomials (`degree = 2` by
#' default), `q > 5` is required to ensure that the local
#' smoothing kernel assigns at least 3 observations a positive
#' weight. This will guarantee a unique least squares fit.
#' When fitting local linear polynomials (`degree = 1`),
#' `q > 3` is sufficient.
#'
#' To avoid k-d tree warnings, set the optional argument
#' `control = loess.control(surface = "direct")`.
#' See [stats::loess.control()] and
#' \insertCite{ClevGros+91;textual}{fastbeta} for details.
#'
#' @param formula A formula specifying a time series in `data`.
#'   Something like `x ~ t`.
#' @param data A data frame containing time series data.
#'   Must satisfy `all(all.vars(formula) %in% names(data))`
#'   (column names match variable names),
#'   `length(unique(diff(data[[all.vars(formula)[2]]]))) == 1`
#'   (time column is equally spaced).
#' @param q An integer vector. The values of `span` passed to
#'   [stats::loess()] will be `q / nrow(data)` (see Details).
#' @param ... Additional arguments of [stats::loess()], such as
#'   `degree`, but not `subset` or `span`. Those not specified
#'   are assigned their default values in [stats::loess()].
#' @param plot A logical scalar. If `TRUE`, then the fitted curves
#'   are plotted (4 curves per plot).
#'
#' @return
#' Invisibly, a list of loess objects with names `paste0("q", q)`.
#'
#' @examples
#' times <- seq(0, 4, by = 0.01)
#' x <- sin(2 * pi * times) + rnorm(times, 0, 1)
#' df <- data.frame(t = times, x)
#' q <- c(6, 20, 100, 400)
#' control <- loess.control(surface = "direct")
#' try_loess_out <- try_loess(x ~ t, df, q, control = control, plot = TRUE)
#'
#' # Fit with `q = 100` looks okay
#' my_loess <- try_loess_out[["q100"]]
#'
#' @references
#' \insertRef{ClevGros+91}{fastbeta}
#'
#' @export
#' @importFrom stats loess
#' @importFrom graphics par plot.new plot.window axis lines title box
try_loess <- function(formula, data, q = 6:9, ..., plot = TRUE) {
  span <- q / nrow(data)
  out <- lapply(span, function(x) loess(formula, data, span = x, ...))
  names(out) <- paste0("q", q)
  if (plot) {
    op <- par(mfrow = c(2, 2), mar = c(3, 4, 2, 1), mgp = c(3, 0.7, 0))
    on.exit(par(op))
    for (i in seq_along(q)) {
      plot.new()
      plot.window(
        xlim = range(out[[i]]$x),
        ylim = range(out[[i]]$y),
        xaxs = "i"
      )
      axis(side = 1)
      axis(side = 2, las = 1)
      lines(formula, data, lwd = 2, col = "grey80")
      lines(out[[i]]$x, out[[i]]$fitted, lwd = 2, col = "seagreen")
      title(main = paste("q =", q[i]), line = 0.5, cex.main = 0.9)
      box()
    }
  }
  invisible(out)
}
