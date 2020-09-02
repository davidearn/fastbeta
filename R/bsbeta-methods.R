#' Methods for class "bsbeta"
#'
#' Method for plotting bsbeta objects returned by [bsbeta()].
#'
#' @param x A bsbeta object.
#' @param ... Unused optional arguments.
#'
#' @name bsbeta-methods
NULL

#' @rdname bsbeta-methods
#' @export
#' @importFrom graphics par plot.new plot.window axis title polygon lines box
plot.bsbeta <- function(x, ...) {
  if (!inherits(x, "bsbeta")) {
    stop("`x` must be a bsbeta object.")
  }
  op <- par(mar=c(3.7,5.5,0.2,0.2)+1, mgp=c(3,0.7,0))
  on.exit(par(op))
  times <- seq_along(x$beta)-1
  ylim <- range(x$beta, x$ci95, na.rm=TRUE)
  plot.new()
  plot.window(xlim=range(times), ylim=ylim, xaxs="i")
  axis(side=1)
  axis(side=2, las=1)
  title(xlab="Time (units dt)", line=2.5)
  title(ylab="Transmission rate", line=4.5)
  index_no_na <- apply(x$ci95, 1, function(x) all(is.finite(x)))
  polygon(x=c(times[index_no_na],rev(times[index_no_na])),
          y=c(x$ci95[index_no_na, 1], rev(x$ci95[index_no_na, 2])),
          col="grey80",
          border=NA)
  lines(times, x$beta)
  box()
  invisible(NULL)
}
