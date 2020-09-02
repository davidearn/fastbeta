#' Methods for class "deconvol"
#'
#' Method for plotting deconvol objects returned by [deconvol()].
#'
#' @param x A deconvol object.
#' @param ... Unused optional arguments.
#'
#' @name deconvol-methods
NULL

#' @rdname deconvol-methods
#' @export
#' @importFrom graphics par plot.new plot.window lines axis title box legend
plot.deconvol <- function(x, ...) {
  if (!inherits(x, "deconvol")) {
    stop("`x` must be a deconvol object.")
  }
  n <- ncol(x$inc)
  op <- par(mar=c(4.1,2.1,2.1,0.5)+1)
  on.exit(par(op))
  ylim <- range(x$x_pad, x$inc, na.rm=TRUE)
  ylim <- c(max(ylim[1]-0.04*diff(ylim), 0), ylim[2]+0.04*diff(ylim))
  plot.new()
  plot.window(xlim=range(x$times), ylim=ylim, xaxs="i", yaxs="i")
  axis(side=1, mgp=c(3,0.7,0))
  axis(side=2, mgp=c(3,0.7,0), las=1)
  title(xlab="Time (units dt)")
  with(x, {
    for (j in seq_len(n)) {
      lines(times, inc[, j], lwd=2, col="grey80")
    }
    lines(times, inc[, n], lwd=2)
    if (any(isTRUE(chi2 < 1))) {
      index_chi2lt1 <- min(which(chi2 < 1))
      lines(times, inc[, index_chi2lt1], lwd=2, col="hotpink")
    }
    lines(times, x_pad, lwd=2, col = "blue")
  })
  box()
  legend("topleft", bty="n", xpd=NA,
         legend=c("Deconvolved incidence (intermediate)",
                  "Deconvolved incidence (final)",
                  "Reported incidence",
                  "Deconvolved incidence (chi2 < 1)"),
         lwd=2, col=c("grey80","black","blue","hotpink"),
         ncol=2, cex=0.8, inset=c(0,-0.15))

  ylim <- range(x$x_pad, x$inc_rep, na.rm=TRUE)
  ylim <- c(max(ylim[1]-0.04*diff(ylim), 0), ylim[2]+0.04*diff(ylim))
  plot.new()
  plot.window(xlim=range(x$times), ylim=ylim, xaxs="i", yaxs="i")
  axis(side=1, mgp=c(3,0.7,0))
  axis(side=2, mgp=c(3,0.7,0), las=1)
  title(xlab="Time (units dt)")
  with(x, {
    for (j in seq_len(n)) {
      lines(times, inc_rep[, j], lwd=2, col="grey80")
    }
    lines(times, inc_rep[, n], lwd=2)
    if (any(isTRUE(chi2 < 1))) {
      lines(times, inc_rep[, index_chi2lt1], lwd=2, col="hotpink")
    }
    lines(times, x_pad, lwd=2, col = "blue")
  })
  box()
  legend("topleft", bty="n", xpd=NA,
         legend=c("Expected reported incidence (intermediate)",
                  "Expected reported incidence (final)",
                  "Reported incidence",
                  "Expected reported incidence (chi2 < 1)"),
         lwd=2, col=c("grey80","black","blue","hotpink"),
         ncol=2, cex=0.8, inset=c(0,-0.15))

  par(mar=c(4,4,0.2,0.2)+1)
  plot(0:(n-1), x$chi2, las=1, xlab="Iterations", ylab="chi2")
  abline(h = 1, lty=2)
  invisible(NULL)
}
