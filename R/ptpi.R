#' \loadmathjax
#' Estimate \mjseqn{S_0} using PTPI
#'
#' @description
#' Using peak-to-peak iteration (PTPI, see Algorithm), estimates the
#' initial number of susceptible individuals \mjseqn{S_0 = S(t_0)}
#' from time series of incidence, births, and natural mortality with
#' observations at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}.
#' The incidence time series must be roughly periodic.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#' Supplied in the data frame `df` are time series
#' \mjseqn{\lbrace (t_i,Z_i) \rbrace_{i=0}^{n-1}},
#' \mjseqn{\lbrace (t_i,B_i) \rbrace_{i=0}^{n-1}}, and
#' \mjseqn{\lbrace (t_i,\mu_i) \rbrace_{i=0}^{n-1}}
#' of incidence, births, and natural mortality with observations
#' at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}.
#' The first peak in incidence and the last peak in phase with
#' the first peak occur at times \mjseqn{t_a} and \mjseqn{t_b},
#' respectively, where \mjseqn{a} is `peak1 - 1` and \mjseqn{b}
#' is `peak2 - 1`.
#'
#' Let \mjseqn{S_i^{(j)}} denote the estimate of \mjseqn{S(t_i)}
#' obtained after \mjseqn{j} iterations. Assign both \mjseqn{S_0^{(0)}}
#' and \mjseqn{S_a^{(0)}} the value of `S0_init`, then proceed with
#' the following iteration.
#'
#' Given \mjseqn{S_a^{(j)}} for some \mjseqn{j}, define \mjseqn{S_i^{(j)}}
#' for all \mjseqn{i \in \lbrace a+1,\ldots,n-1 \rbrace} with the recursion
#'
#' \mjsdeqn{S_{i+1}^{(j)} = \frac{\big\lbrack 1 - \frac{1}{2} \mu_i \Delta t \big\rbrack S_i^{(j)} + B_{i+1} - Z_{i+1}}{1 + \frac{1}{2} \mu_{i+1} \Delta t}\,.}
#'
#' Then set \mjseqn{S_a^{(j+1)} = S_b^{(j)}} and define \mjseqn{S_i^{(j+1)}}
#' for all \mjseqn{i \in \lbrace 0,\ldots,a-1 \rbrace} with the recursion
#'
#' \mjsdeqn{S_{i-1}^{(j+1)} = \frac{\big\lbrack 1 + \frac{1}{2} \mu_i \Delta t \big\rbrack S_i^{(j+1)} - B_i + Z_i}{1 + \frac{1}{2} \mu_{i-1} \Delta t}\,.}
#'
#' Note that this second recursion is carried out backwards in time.
#'
#' ## 2. Convergence
#' Suppose incidence is \mjseqn{(p \Delta t)}-periodic, where \mjseqn{p}
#' is a nonzero integer, and consider the special case in which the birth
#' and per capita death rates are constant. Then \mjseqn{Z_i = Z_{i+kp}},
#' \mjseqn{B_i = B}, and \mjseqn{\mu_i = \mu}, for all integers \mjseqn{i}
#' and \mjseqn{k}. In this case, one can show that the sequence
#' \mjseqn{(S_0^{(j)})} converges to a limit \mjseqn{S_0^*} given by
#'
#' \mjsdeqn{S_0^* = \left( \frac{1 + \frac{1}{2} \mu \Delta t}{1 - \frac{1}{2} \mu \Delta t} \right)^a S_a^* + \frac{1}{1 - \frac{1}{2} \mu \Delta t} \sum_{i=1}^a (Z_i - B_i) \left( \frac{1 + \frac{1}{2} \mu \Delta t}{1 - \frac{1}{2} \mu \Delta t} \right)^{i-1}\,,}
#'
#' where \mjseqn{S_a^*} is the limit of the sequence \mjseqn{(S_a^{(j)})},
#' given by
#'
#' \mjsdeqn{S_a^* = \frac{B}{\mu \Delta t} - \frac{1}{1 + \frac{1}{2} \mu \Delta t} \left\lbrack \frac{\left( \frac{1 - \frac{1}{2} \mu \Delta t}{1 + \frac{1}{2} \mu \Delta t} \right)^p}{1 - \left( \frac{1 - \frac{1}{2} \mu \Delta t}{1 + \frac{1}{2} \mu \Delta t} \right)^p} \right\rbrack \sum_{i=1}^p Z_{a+i} \left( \frac{1 + \frac{1}{2} \mu \Delta t}{1 - \frac{1}{2} \mu \Delta t} \right)^i\,.}
#'
#' One can further show that the convergence is linear and follows
#'
#' \mjsdeqn{S_0^{(j)} - S_0^* = (S_0^{(0)} - S_a^*) \left( \frac{1 - \frac{1}{2} \mu \Delta t}{1 + \frac{1}{2} \mu \Delta t} \right)^{jp-a}\,.}
#'
#' @param df A data frame with numeric columns:
#'
#'   \describe{
#'     \item{`Z`}{Incidence. `Z[i]` is the number of infections between
#'       times `t[i-1]` and `t[i]`. Must be roughly periodic.
#'     }
#'     \item{`B`}{Births. `B[i]` is the number of births between times
#'       `t[i-1]` and `t[i]`.
#'     }
#'     \item{`mu`}{Natural mortality rate. `mu[i]` is the rate at time
#'       `t[i]` expressed per unit \mjseqn{\Delta t} and per capita.
#'     }
#'   }
#'
#'   Missing values in `df` are not tolerated and must be imputed
#'   separately.
#' @param peak1 Integer scalar. Index of the first peak in `df$Z`.
#'   A reasonable value can be obtained using [peaks()] (see Examples).
#' @param peak2 Integer scalar. Index of the last peak in `df$Z` in phase
#'   with the first peak.
#'   A reasonable value can be obtained using [peaks()] (see Examples).
#' @param S0_init Numeric scalar. An initial estimate of \mjseqn{S_0}.
#' @param it Integer scalar. The number of iterations to perform before
#'   stopping.
#'
#' @return
#' A ptpi object. A list with elements:
#'
#' \describe{
#'   \item{`mat`}{A numeric matrix with `nrow(df)` rows and `1 + it` columns
#'     containing the susceptible time series generated in each iteration.
#'     `mat[i,j]` gives the value of \mjseqn{S_{i-1}^{(j-1)}} as defined in
#'     Algorithm.
#'   }
#'   \item{`S0`}{A numeric vector listing in order all `1 + it` estimates of
#'     \mjseqn{S_0 = S(t_0)}. Equivalent to `mat[1, ]`.
#'   }
#'   \item{`S0_final`}{The final estimate of \mjseqn{S_0}.
#'     Equivalent to `mat[1, ncol(mat)]`.
#'   }
#' }
#'
#' The object has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(ptpi, arg_list)`.
#'
#' @examples
#' # Simulate 20 years of disease incidence,
#' # observed weekly
#' par_list <- make_par_list(dt_weeks = 1)
#' df <- make_data(
#'   par_list = par_list,
#'   with_ds = TRUE,
#'   n = 20 * 365 / 7 # number of weeks in 20 years
#' )
#'
#' # Plot incidence time series, and note the
#' # apparent 1-year period
#' plot(Z ~ I(t_years - t_years[1]), df, type = "l",
#'   xlab = "Time (years)",
#'   ylab = "Incidence"
#' )
#'
#' # Locate peaks in incidence time series
#' peaks_out <- peaks(
#'   x = df$Z,
#'   bw1 = 6,
#'   bw2 = 8,
#'   period = 365 / 7 # number of weeks in 1 year
#' )
#'
#' # Verify that peaks were identified correctly:
#' # * Blue lines are peaks.
#' # * Pink circles are peaks in phase with first peak.
#' plot(peaks_out)
#'
#' # Estimate `S0` from incidence, births, and
#' # natural mortality using PTPI, starting from
#' # an erroneous initial guess
#' ptpi_out <- ptpi(
#'   df = df,
#'   S0_init = df$S[1] * 4,
#'   peak1 = with(peaks_out, phase[1]),
#'   peak2 = with(peaks_out, phase[length(phase)]),
#'   it = 25
#' )
#' plot(ptpi_out)
#'
#' # Sequence of estimates
#' ptpi_out$S0
#'
#' # Relative error in final estimate
#' (ptpi_out$S0_final - df$S[1]) / df$S[1]
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' @seealso [peaks()], [plot.ptpi()]
#' @export
ptpi <- function(df, S0_init, peak1 = 1L, peak2 = nrow(df), it = 10L) {
  # Save arguments in a list
  arg_list <- as.list(environment())

  # Preallocate memory for all susceptible time series,
  # and initialize the first
  mat <- matrix(NA, nrow = nrow(df), ncol = 1 + it)
  mat[1, 1] <- S0_init
  mat[peak1, 1] <- S0_init

  # Iterate
  for (j in seq_len(1 + it)) {

    # Update estimate of susceptibles at index `peak1`
    # with estimate at index `peak2` after reconstructing
    # from index `peak1` to end
    for (i in (peak1+1):nrow(mat)) {
      mat[i, j] <- with(df[c("Z", "B", "mu")],
        {
          ((1 - 0.5 * mu[i-1] * 1) * mat[i-1,j] + B[i] - Z[i]) /
            (1 + 0.5 * mu[i] * 1)
        }
      )
    }
    if (j == it + 1) {
      break
    }
    mat[peak1, j+1] <- mat[peak2, j]

    # Update estimate of susceptibles at initial time
    # by reconstructing from index `peak1` to start
    # (backwards in time)
    for (i in (peak1-1):1) {
      mat[i, j+1] <- with(df[c("Z", "B", "mu")],
        {
          ((1 + 0.5 * mu[i+1] * 1) * mat[i+1, j+1] - B[i+1] + Z[i+1]) /
            (1 - 0.5 * mu[i] * 1)
        }
      )
    }

  }

  out <- list(
    mat      = mat,
    S0       = mat[1, ],
    S0_final = mat[1, ncol(mat)]
  )
  structure(out,
    class    = c("ptpi", "list"),
    call     = match.call(),
    arg_list = arg_list
  )
}

#' Methods for class ptpi
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
       xaxs="i", las=1,
       xlab="Time (units dt)", ylab="")
  mtext("Susceptibles", side=2, line=4.5)
  for (j in seq_len(n-1)) {
    lines(0:(m-1), x$mat[, j], lwd=2, col="grey80")
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
