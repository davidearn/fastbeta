#' Methods for class "fastbeta"
#'
#' Methods for printing and plotting fastbeta objects
#' returned by [fastbeta()].
#'
#' @param x A fastbeta object.
#' @param ... Unused optional arguments.
#'
#' @name fastbeta-methods
NULL

#' @rdname fastbeta-methods
#' @export
#' @import graphics
plot.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  ynames <- c("S", "I", "beta")
  ylabs <- c("Susceptible", "Infectious", "Transmission rate")
  cols <- c("seagreen", "mediumvioletred", "slateblue")
  op <- par(mfrow=c(3,1), mar=c(3,5.6,0.2,0.2), oma=c(1.1,0,1.9,0)+1)
  on.exit(par(op))
  for (i in 1:3) {
    plot(seq_len(nrow(x$out))-1, x$out[, ynames[i]],
         type="l", lty=1, lwd=2, col=cols[i],
         xaxs="i", las=1, mgp=c(3,0.7,0),
         xlab="", ylab="")
    title(ylab=ylabs[i], line=4.5, cex.lab=1.2)
    if (i == 1) {
      title(main=paste(toupper(x$method), "method"), line=1, cex.main=1.5, xpd=NA)
    }
    if (i == 3) {
      title(xlab="Time (units dt)", line=3, cex.lab=1.5, xpd=NA)
    }
  }
  invisible(NULL)
}

#' @rdname fastbeta-methods
#' @export
print.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  n <- nrow(x$out)
  print(x$out[seq_len(min(n, 16)), ])
  if (n > 16) {
    s <- if (n == 17) "" else "s"
    cat("\n... and ", n - 16, " more row", s, ".\n", sep = "")
  }
  invisible(x$out)
}
