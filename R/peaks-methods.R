#' Methods for class "peaks"
#'
#' Methods for plotting and printing peaks objects
#' returned by [peaks()].
#'
#' @param x A peaks object.
#' @param ... Unused optional arguments.
#'
#' @name peaks-methods
NULL

#' @rdname peaks-methods
#' @export
#' @import graphics
plot.peaks <- function(x, ...) {
  if (!inherits(x, "peaks")) {
    stop("`x` must be a peaks object.")
  }
  op <- par(mar=c(4,4.4,0.8,1)+1)
  on.exit(par(op))
  plot(seq_along(x$x)-1, x$x, type="l",
       lwd=2, col="grey70",
       xaxs="i", las=1, mgp=c(3,0.7,0),
       xlab="Time (units dt)", ylab="")
  mtext("x", side=2, line=4, las=1)
  lines(seq_along(x$x)-1, x$xbar, lwd=2)
  abline(v=x$all-1, col="blue")
  if (!is.null(x$phase)) {
    axis(side=3, at=x$phase-1,
         labels=rep("o", length(x$phase)),
         tick=FALSE, mgp=c(3,0.1,0), col.axis="hotpink")
  }
  invisible(NULL)
}

#' @rdname peaks-methods
#' @export
print.peaks <- function(x, ...) {
  if (!inherits(x, "peaks")) {
    stop("`x` must be a peaks object.")
  }
  n <- length(x$all)
  if (n == 0) {
    cat("`peaks()` found no peaks.\n")
  } else {
    fmt <- paste0("%", nchar(max(x$all)) + 4, "d")
    cat("`peaks()` found", n, "peaks in `x`, indexed by:\n")
    cat(paste0(sprintf(fmt, x$all), "\n"), sep = "")
  }
  invisible(x$all)
}
