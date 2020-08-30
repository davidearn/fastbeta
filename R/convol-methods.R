#' Methods for class "convol"
#'
#' Method for plotting convol objects returned by [convol()].
#'
#' @param x A convol object.
#' @param ... Unused optional arguments.
#'
#' @name convol-methods
NULL

#' @rdname convol-methods
#' @export
#' @importFrom graphics par plot.new plot.window lines axis title box legend
plot.convol <- function(x, ...) {
  if (!inherits(x, "convol")) {
    stop("`x` must be a convol object.")
  }
  op <- par(mar=c(4.1,2.1,2.1,0.5)+1)
  on.exit(par(op))
  ylim <- range(x$x_pad, x$convolution, x$simulation, na.rm=TRUE)
  ylim <- c(max(ylim[1]-0.04*diff(ylim), 0), ylim[2]+0.04*diff(ylim))
  plot.new()
  plot.window(xlim=range(x$times), ylim=ylim, xaxs="i", yaxs="i")
  axis(side=1, mgp=c(3,0.7,0))
  axis(side=2, mgp=c(3,0.7,0), las=1)
  title(xlab="Time (units dt)")
  with(x, {
    if (!is.null(simulation)) {
      n <- ncol(simulation)
      for (j in seq_len(n)) {
        lines(times, simulation[, j], lwd=2, col="grey80")
      }
    }
    lines(times, convolution, lwd=2)
    lines(times, x_pad, lwd=2, col="blue")
  })
  box()
  legend("topleft", bty="n", xpd=NA,
         legend=c("Reported incidence (simulated)",
                  "Reported incidence (expected)",
                  "Incidence"),
         lwd=2, col=c("grey80","black","blue"),
         ncol=2, cex=0.8, inset=c(0,-0.15))
  invisible(NULL)
}
