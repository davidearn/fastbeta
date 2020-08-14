#' Plot a fastbeta object
#'
#' @description
#' Plots a fastbeta object.
#'
#' @param x A fastbeta object.
#' @param ... Unused optional arguments.
#'
#' @return
#' `NULL` is returned invisibly.
#'
#' @export
plot.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  graphics::layout(matrix(1:3, ncol = 1))
  ynames <- c("S", "I", "beta")
  ylabs <- c(
    "Number of susceptibles",
    "Number of infecteds",
    "Transmission rate, per unit dt\nper susceptible per infected"
  )
  cols <- c("green", "red", "blue")
  for (i in 1:3) {
    graphics::plot(0:(nrow(x)-1), x[, ynames[i]],
      type = "l",
      lty  = 1,
      lwd  = 2,
      col  = cols[i],
      xaxs = "i",
      las  = 1,
      xlab = "Time, units dt",
      ylab = ylabs[i]
    )
  }
  invisible(NULL)
}
