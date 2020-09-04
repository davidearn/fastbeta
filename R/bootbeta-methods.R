#' Methods for class "bootbeta"
#'
#' Method for plotting bootbeta objects returned by [bootbeta()].
#'
#' @param x A bootbeta object.
#' @param ... Unused optional arguments.
#'
#' @name bootbeta-methods
NULL

#' @rdname bootbeta-methods
#' @export
#' @import graphics
plot.bootbeta <- function(x, ...) {
  if (!inherits(x, "bootbeta")) {
    stop("`x` must be a bootbeta object.")
  }
  op <- par(mar=c(3.7,5.5,0.2,0.2)+1, mgp=c(3,0.7,0))
  on.exit(par(op))
  ylim <- range(x$beta, x$ci95, na.rm=TRUE)
  plot.new()
  plot.window(xlim=range(x$times), ylim=ylim, xaxs="i")
  axis(side=1)
  axis(side=2, las=1)
  title(xlab="Time (units dt)", line=2.5)
  title(ylab="Transmission rate", line=4.5)
  index_no_na <- apply(x$ci95, 1, function(x) all(is.finite(x)))
  polygon(x=c(x$times[index_no_na],rev(x$times[index_no_na])),
          y=c(x$ci95[index_no_na, 1], rev(x$ci95[index_no_na, 2])),
          col="grey80",
          border=NA)
  lines(x$times, x$beta)
  box()
  invisible(NULL)
}
