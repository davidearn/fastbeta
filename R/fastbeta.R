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
#' @param method A character scalar, one of `"fc"`, `"s"`, `"si"`,
#'   and `"sei"`, indicating a transmission rate estimation method.
#'   The FC and S methods are deprecated and outperformed by the
#'   more robust SI and SEI methods. For the details of each algorithm,
#'   see the respective sections of [this help page][estimate-beta].
#'
#' @return
#' A fastbeta object.
#'
#' @examples
#' # Simulate time series data using an SIR model
#' # with seasonally forced transmission rate
#' par_list <- make_par_list(model = "sir")
#' df <- make_data(par_list = par_list, with_ds = TRUE, model = "sir")
#' 
#' # Estimate the seasonally forced transmission rate
#' # using the SI method
#' par_list <- within(par_list, tgen <- tlat + tinf)
#' fastbeta_out <- fastbeta(df, par_list, method = "si")
#' plot(fastbeta_out)
#' 
#' # Simulate time series data using an SEIR model
#' # with seasonally forced transmission rate
#' par_list <- make_par_list(model = "seir")
#' df <- make_data(par_list = par_list, with_ds = TRUE, model = "seir")
#' 
#' # Estimate the seasonally forced transmission rate
#' # using the SEI method
#' fastbeta_out <- fastbeta(df, par_list, method = "sei")
#' plot(fastbeta_out)
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' @seealso [methods for class "fastbeta"][fastbeta-methods],
#'   [algorithms used under the hood][estimate-beta]
#' @export
fastbeta <- function(df, par_list, method) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  } else if (!all(c("t", "Z", "B", "mu") %in% names(df))) {
    stop("`df` is missing necessary columns.")
  }
  if (!is.list(par_list)) {
    stop("`par_list` must be a list.")
  } else if ((method %in% c("fc", "s") &&
                !all(c("S0", "tgen") %in% names(par_list))) ||
             (method == "si" &&
                !all(c("S0", "I0", "tgen") %in% names(par_list))) ||
             (method == "sei" &&
                !all(c("S0", "E0", "I0", "tlat", "tinf") %in% names(par_list)))) {
    stop("`par_list` is missing necessary elements.")
  }

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Missing values are not tolerated. Zeros in incidence
  # *are* tolerated but can have undesired effects. See
  # `?"estimate-beta"` for details.
  df[c("Z", "B", "mu")] <- mapply(impute_na,
    x = df[c("Z", "B", "mu")],
    zero_as_na = c(TRUE, FALSE, FALSE)
  )

  # Apply the desired method
  f <- get(paste0("estimate_beta_", method))
  out <- f(df, par_list)

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
#' @importFrom graphics layout par plot title mtext
plot.fastbeta <- function(x, ...) {
  if (!inherits(x, "fastbeta")) {
    stop("`x` must be a fastbeta object.")
  }
  method <- attr(x, "arg_list")$method
  layout(matrix(1:3, ncol = 1))
  ynames <- c("S", "I", "beta")
  ylabs <- c(
    "Susceptible",
    if (method == "sei") "Infectious" else "Infected",
    "Transmission rate"
  )
  cols <- c("seagreen", "mediumvioletred", "slateblue")
  op <- par(mar=c(3,5.6,0.2,0.2), oma=c(1.1,0,1.9,0)+1)
  on.exit(par(op))
  for (i in 1:3) {
    plot(seq_len(nrow(x))-1, x[, ynames[i]],
         type="l", lty=1, lwd=2, col=cols[i],
         xaxs="i", las=1, xlab="", ylab="")
    title(ylab=ylabs[i], line=4.5, cex.lab=1.2)
    if (i == 1) {
      title(main=paste(toupper(method), "method"), line=1, cex.main=1.5, xpd=NA)
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
  n <- nrow(x)
  print.data.frame(x[seq_len(min(n, 16)), ])
  if (n > 16) {
    s <- if (n == 17) "" else "s"
    cat("\n... and ", n - 16, " more row", s, ".\n", sep = "")
  }
  invisible(x)
}
