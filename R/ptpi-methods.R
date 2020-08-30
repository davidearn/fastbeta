#' Methods for class "ptpi"
#'
#' Methods for plotting and printing ptpi objects
#' returned by [ptpi()].
#'
#' @param x A ptpi object.
#' @param ... Unused optional arguments.
#'
#' @name ptpi-methods
NULL

#' @rdname ptpi-methods
#' @export
#' @importFrom graphics plot lines points mtext
plot.ptpi <- function(x, ...) {
  if (!inherits(x, "ptpi")) {
    stop("`x` must be a ptpi object.")
  }
  m <- nrow(x$mat)
  n <- ncol(x$mat)
  if (m < 2 || n < 1) {
    stop("`x$mat` must have at least 2 rows and at least 1 column.")
  }
  op <- par(mar=c(4,5.4,1.2,0.2)+1)
  on.exit(op)
  plot(0, 0, type="n",
       xlim=c(0,m-1), ylim=range(x$mat, na.rm=TRUE),
       xaxs="i", las=1, mgp=c(3,0.7,0),
       xlab="Time (units dt)", ylab="")
  mtext("Susceptibles", side=2, line=4.5)
  for (j in seq_len(n-1)) {
    lines(0:(m-1), x$mat[, j], lwd=2, col="grey70")
  }
  lines(0:(m-1), x$mat[, n], lwd=2)
  points(rep(0, length(x$S0)), x$S0, pch=18, col="seagreen", xpd=NA)
  mtext(paste("S0 =", sprintf("%.4f", x$S0_final), "after", n-1, "iterations"),
        side=3, line=0.2, adj=0, padj=0)
  invisible(NULL)
}

#' @rdname ptpi-methods
#' @export
print.ptpi <- function(x, ...) {
  if (!inherits(x, "ptpi")) {
    stop("`x` must be a ptpi object.")
  }
  n <- length(x$S0)
  if (n == 1) {
    cat("`ptpi()` returned this 1 estimate of S0:\n")
  } else {
    cat("`ptpi()` returned these", length(x$S0), "estimates of S0:\n")
  }
  fmt <- paste0("%", nchar(max(floor(x$S0))) + 4 + 5, ".4f")
  cat(paste0(sprintf(fmt, x$S0), "\n"), sep = "")
  invisible(x$S0)
}
