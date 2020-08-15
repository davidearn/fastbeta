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
#' @importFrom graphics layout par plot title mtext
#' @export
plot.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  layout(matrix(1:3, ncol = 1))
  ynames <- c("S", "I", "beta")
  ylabs <- c("Susceptibles", "Infecteds", "Transmission rate")
  cols <- c("seagreen", "mediumvioletred", "slateblue")
  op <- par(mar=c(3, 6, 1, 1), oma=c(3, 0, 0, 0))
  on.exit(par(op))
  for (i in 1:3) {
    plot(0:(nrow(x)-1), x[, ynames[i]],
         type="l", lty=1, lwd=2, col=cols[i],
         xaxs="i", las=1, xlab="", ylab="")
    title(ylab=ylabs[i], line=4.5, cex.lab=1.3)
    mtext("Time (units dt)", side=1, line=1.5, adj=0.565, outer=TRUE, cex=1.3)
  }
  invisible(NULL)
}
