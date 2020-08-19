#' Estimate time-varying transmission rates
#'
#' @description
#' Generates a fastbeta object.
#'
#' @details
#' Details to come.
#'
#' @param df A data frame with time series data.
#' @param par_list A list of parameter values.
#'
#' @return
#' A fastbeta object.
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#' 
#' @seealso [estimate_beta_si()], [deconvolve()], [plot.fastbeta()] 
#' @export
fastbeta <- function(df, par_list) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  } else if (!all(c("t", "Z", "B", "mu") %in% names(df))) {
    stop("`df` is missing necessary columns.")
  }
  if (!is.list(par_list)) {
    stop("`par_list` must be a list.")
  } else if (!all(c("S0", "I0", "tgen") %in% names(par_list))) {
    stop("`par_list` is missing necessary elements.")
  }

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Missing values are not tolerated. Zeros in incidence
  # *are* tolerated but can have undesired effects. See
  # `?estimate_beta_si`.
  df[c("Z", "B", "mu")] <- mapply(impute_na,
    x = df[c("Z", "B", "mu")],
    zero_as_na = c(TRUE, FALSE, FALSE)
  )
  out <- estimate_beta_si(df, par_list)

  # Negative susceptibles indicates that incidence was
  # overestimated or births were underestimated or both
  if (any(out$S < 0, na.rm = TRUE)) {
    warning(
      "Negative elements in susceptibles (`S`) column. Retry with:",
      "\n* scaled down incidence (`df$Z`), and/or",
      "\n* scaled up births (`df$B`)."
    )
  }

  structure(out,
    class = c("fastbeta", "data.frame"),
    call = match.call(),
    arg_list = arg_list
  )
}

#' Methods for class fastbeta
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
#' @importFrom graphics layout par plot title mtext
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
    plot(seq_len(nrow(x))-1, x[, ynames[i]],
         type="l", lty=1, lwd=2, col=cols[i],
         xaxs="i", las=1, xlab="", ylab="")
    title(ylab=ylabs[i], line=4.5, cex.lab=1.3)
    mtext("Time (units dt)", side=1, line=1.5, adj=0.565, outer=TRUE, cex=1.3)
  }
  invisible(NULL)
}

#' @rdname fastbeta-methods
#' @export
print.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  n <- nrow(x)
  print.data.frame(x[seq_len(min(n, 16)), ])
  if (n > 16) {
    s <- if (n == 17) "" else "s"
    cat("\n... and ", n - 16, " more row", s, ".\n", sep = "")
  }
  invisible(x)
}
