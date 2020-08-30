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
#' Supplied in the data frame `data` are time series
#' \mjseqn{\lbrace (t_i,Z_i) \rbrace_{i=0}^{n-1}},
#' \mjseqn{\lbrace (t_i,B_i) \rbrace_{i=0}^{n-1}}, and
#' \mjseqn{\lbrace (t_i,\mu_i) \rbrace_{i=0}^{n-1}}
#' of incidence, births, and natural mortality with observations
#' at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}.
#' The first peak in incidence and the last peak in phase with
#' the first peak occur at times \mjseqn{t_a} and \mjseqn{t_b},
#' respectively, where \mjseqn{a} is `peak1-1` and \mjseqn{b}
#' is `peak2-1`.
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
#' Suppose that incidence is \mjseqn{(p \Delta t)}-periodic and
#' (without loss of generality) that \mjseqn{b - a = m p}, where
#' \mjseqn{p} and \mjseqn{m} are positive integers. Then
#' \mjseqn{Z_i = Z_{i+kp}} for all integers \mjseqn{i} and \mjseqn{k},
#' and there are \mjseqn{m} periods between times \mjseqn{t_a} and
#' \mjseqn{t_b}. Now consider the special case in which the birth
#' and per capita death rates are constant, so that \mjseqn{B_i = B}
#' and \mjseqn{\mu_i = \mu} for all integers \mjseqn{i}.
#' In this case, one can show that the sequence \mjseqn{(S_0^{(j)})}
#' converges to a limit \mjseqn{S_0^*} given by
#'
#' \mjsdeqn{S_0^* = \left( \frac{1 + \frac{1}{2} \mu \Delta t}{1 - \frac{1}{2} \mu \Delta t} \right)^a S_a^* + \frac{1}{1 - \frac{1}{2} \mu \Delta t} \sum_{i=1}^a (Z_i - B) \left( \frac{1 + \frac{1}{2} \mu \Delta t}{1 - \frac{1}{2} \mu \Delta t} \right)^{i-1}\,,}
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
#' A corollary is that convergence of \mjseqn{(S_0^{(j)})} is guaranteed
#' even if incidence is not truly periodic, because, within the iteration,
#' incidence is always *effectively* \mjseqn{((b - a) \Delta t)}-periodic.
#' This can be seen by noting that \mjseqn{S_a^{(j)}} is updated recursively,
#' with each update using the same \mjseqn{b - a} observations of incidence:
#' \mjseqn{Z_{a+1}, \ldots, Z_b}. Hence the limits above will be correct for
#' aperiodic (or just roughly periodic) incidence provided that \mjseqn{p}
#' is replaced with \mjseqn{b - a}.
#'
#' @param data A data frame with numeric columns:
#'
#'   \describe{
#'     \item{`t`}{\mjseqn{\lbrace\,t_i\,\rbrace}
#'       Time. Lists the observation times \mjseqn{t_i = t_0 + i \Delta t}
#'       (for \mjseqn{i = 0, \ldots, n-1}) in any convenient units. Here,
#'       \mjseqn{t_0} is the initial observation time and \mjseqn{\Delta t}
#'       is the observation interval.
#'     }
#'     \item{`Z`}{\mjseqn{\lbrace\,Z_i\,\rbrace}
#'       Incidence. `Z[i]` is the number of infections
#'       between times `t[i-1]` and `t[i]`. Must be roughly periodic.
#'     }
#'     \item{`B`}{\mjseqn{\lbrace\,B_i\,\rbrace}
#'       Births. `B[i]` is the number of births
#'       between times `t[i-1]` and `t[i]`.
#'     }
#'     \item{`mu`}{\mjseqn{\lbrace\,\mu_i \Delta t\,\rbrace}
#'       Per capita natural mortality rate. `mu[i]` is the rate
#'       at time `t[i]` expressed per unit \mjseqn{\Delta t}.
#'     }
#'   }
#'
#' @param peak1 An integer scalar. Index of the first peak in `data$Z`.
#'   A reasonable value can be obtained using [peaks()] (see Examples).
#' @param peak2 An integer scalar. Index of the last peak in `data$Z`
#'   in phase with the first peak.
#'   A reasonable value can be obtained using [peaks()] (see Examples).
#' @param S0_init A numeric scalar. An initial estimate of \mjseqn{S_0}.
#' @param it An integer scalar. The number of iterations to perform
#'   before stopping.
#'
#' @return
#' A ptpi object. A list with elements:
#'
#' \describe{
#'   \item{`times`}{A numeric vector. The rows of `mat` correspond
#'     to these time points. Identical to `data$t`.
#'   }
#'   \item{`mat`}{A numeric matrix with `length(times)` rows and `1 + it`
#'     columns containing the susceptible time series generated in each
#'     iteration. `mat[i,j]` gives the value of \mjseqn{S_{i-1}^{(j-1)}}
#'     as defined in Algorithm.
#'   }
#'   \item{`S0`}{A numeric vector listing in order all `1 + it` estimates of
#'     \mjseqn{S_0 = S(t_0)}. Equivalent to `mat[1, ]`.
#'   }
#'   \item{`S0_final`}{A numeric scalar giving the final estimate
#'     of \mjseqn{S_0}. Equivalent to `mat[1, ncol(mat)]`.
#'   }
#'   \item{`call`}{The function call. The ptpi object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call.
#'     The ptpi object is reproducible with `do.call(ptpi, arg_list)`.
#'   }
#' }
#'
#' @examples
#' # Simulate 20 years of disease incidence,
#' # observed weekly
#' pl <- make_par_list(model = "sir")
#' df <- make_data(pl,
#'   n = 20 * 365 / 7, # number of weeks in 20 years
#'   with_ds = TRUE,
#'   model = "sir"
#' )
#'
#' # Plot incidence time series, and note the
#' # apparent 1-year period
#' plot(Z ~ t_years, df, type = "l",
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
#'   data = df,
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
#' @seealso [methods for class "ptpi"][ptpi-methods], [peaks()]
#' @export
ptpi <- function(data, S0_init, peak1 = 1L, peak2 = nrow(data), it = 10L) {
  if (missing(data)) {
    stop("Missing argument `data`.")
  } else if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  } else if (!all(c("t", "Z", "B", "mu") %in% names(data))) {
    stop("`data` is missing necessary columns.")
  }
  if (missing(S0_init)) {
    stop("Missing argument `S0_init`.")
  } else if (!is.numeric(S0_init) || length(S0_init) != 1 || !isTRUE(S0_init >= 0)) {
    stop("`S0_init` must be a non-negative numeric scalar.")
  }
  if (missing(peak1)) {
    stop("Missing argument `peak1`.")
  } else if (missing(peak2)) {
    stop("Missing argument `peak2`.")
  } else if (!is.numeric(peak1) || length(peak1) != 1 || !peak1 %in% seq_len(nrow(data)) ||
               !is.numeric(peak2) || length(peak2) != 1 || !peak2 %in% seq_len(nrow(data))) {
    stop("`peak1` and `peak2` must be integer scalars in `seq_len(nrow(data))`.")
  } else if (peak2 < peak1) {
    stop("`peak2` must be greater than `peak1`.")
  }
  if (missing(it)) {
    stop("Missing argument `it`.")
  } else if (!is.numeric(it) || length(it) != 1 || !isTRUE(it >= 0)) {
    stop("`it` must be a non-negative integer scalar.")
  }

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Preallocate memory for all susceptible time series,
  # and initialize the first
  it <- floor(it)
  mat <- matrix(NA, nrow = nrow(data), ncol = 1 + it)
  mat[1, 1] <- S0_init
  mat[peak1, 1] <- S0_init

  # Iterate
  for (j in seq_len(1 + it)) {

    # Update estimate of susceptibles at index `peak1`
    # with estimate at index `peak2` after reconstructing
    # from index `peak1` to end
    for (i in (peak1+1):nrow(mat)) {
      mat[i, j] <- with(data[c("Z", "B", "mu")],
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
      mat[i, j+1] <- with(data[c("Z", "B", "mu")],
        {
          ((1 + 0.5 * mu[i+1] * 1) * mat[i+1, j+1] - B[i+1] + Z[i+1]) /
            (1 - 0.5 * mu[i] * 1)
        }
      )
    }

  }

  out <- list(
    times    = data$t,
    mat      = mat,
    S0       = mat[1, ],
    S0_final = mat[1, ncol(mat)],
    call     = match.call(),
    arg_list = arg_list
  )
  structure(out, class = c("ptpi", "list"))
}
