#' \loadmathjax
#' Estimate time-varying transmission rates
#'
#' @description
#' A wrapper for the [`estimate_beta_method()`][estimate-beta]
#' functions. These implement the FC, S, SI, and SEI methods for
#' estimating time-varying transmission rates \mjseqn{\beta(t)}
#' from time series of incidence, births, and natural mortality
#' observed at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}.
#' Note that the algorithms require incidence, *not* reported
#' incidence. See [deconvol()] for a method of reconstructing
#' incidence from reported incidence.
#'
#' @details
#' # Details
#'
#' ## 1. Models and algorithms
#'
#' See [here][estimate-beta].
#'
#' ## 2. Missing data and zeros in incidence
#'
#' Missing values in `data` are not tolerated by any of the estimation
#' methods. Columns `S` and `beta` in the output will be filled with `NA`
#' after the index of the first encountered missing value.
#'
#' In the FC and S methods, zeros in `data$Z` can cause divide-by-zero
#' errors. In this case, column `beta` in the output may contain `NaN`
#' and `Inf` where an estimate would otherwise be expected. In the SI
#' and SEI methods, strings of zeros in `data$Z` can lead to spurious
#' zeros and large numeric elements in column `beta` of the output.
#'
#' `fastbeta()` applies [impute_na()] to the columns of `data` in an
#' attempt to impute `NA` in `data` and replace zeros in `data$Z` with
#' positive values. [impute_na()] fails for missing values and zeros
#' at the edges of a vector. If a different imputation scheme is
#' desired ([impute_na()] uses linear interpolation), then it must be
#' applied prior to calling `fastbeta()`.
#'
#' ## 3. Negative susceptibles
#'
#' If true incidence is systematically overestimated by `data$Z`
#' or true births are systematically underestimated by `data$B`,
#' then the estimated susceptible population size can grossly
#' underestimate the true number of susceptibles and can even
#' become negative. In the latter case, column `S` in the output
#' will have negative elements.
#'
#' `fastbeta()` warns if the estimated susceptible population size
#' is at any time negative and suggests restarting with scaled down
#' `data$Z` and/or scaled up `data$B`. A rule-of-thumb is to iterate
#' with different scalings until the susceptible time series no
#' longer displays pronounced transient dynamics, while considering
#' whether the correction applied is reasonable.
#'
#' ## 4. Noise in the estimated transmission rate
#'
#' All four estimation methods allow noise in the data (due to
#' process and observation error) to be propagated to transmission
#' rate estimates. To distill temporal  patterns of interest from
#' the noise, the raw transmission rate estimates can be smoothed
#' using, for example, [stats::loess()]. `fastbeta()` does not
#' undertake smoothing, because determining what constitutes a
#' "good" degree of smoothing is not a trivial task and typically
#' involves a graphical inspection of results for different values
#' of a smoothing parameter. [try_loess()] can be used to greatly
#' facilitate the selection of loess fits to noisy transmission
#' rate estimates.
#'
#' ## 5. Confidence intervals on the estimated transmission rate
#'
#' `fastbeta()` does not construct confidence intervals on any
#' of its estimates. The returned fastbeta object can be passed
#' to [bootbeta()] to obtain bootstrap 95% confidence intervals.
#'
#' ## 6. Effective observation interval (FC method only)
#'
#' In the absence of errors, every `round(par_list$tgen)`th row in
#' the FC method output will be complete. The remaining rows will
#' contain `NA` in columns `S`, `I`, `beta`, `Z_agg`, and `B_agg`.
#' This is not a mistake: the rounded generation interval becomes
#' the *effective* observation interval when incidence and births
#' are aggregated over the rounded generation interval (see the
#' FC method algorithm [here][estimate-beta]).
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
#'       between times `t[i-1]` and `t[i]`.
#'     }
#'     \item{`B`}{\mjseqn{\lbrace\,B_i\,\rbrace}
#'       Births. `B[i]` is the number of births
#'       between times `t[i-1]` and `t[i]`.
#'     }
#'     \item{`mu`}{\mjseqn{\lbrace\,\mu_i \Delta t\,\rbrace}
#'       Per capita natural mortality rate. `mu[i]` is the rate
#'       at time `t[i]` expressed per unit \mjseqn{\Delta t}.
#'       S, SI, and SEI methods only.
#'     }
#'   }
#'
#' @param par_list A list with numeric scalar elements:
#'
#'   \describe{
#'     \item{`tgen`}{\mjseqn{\lbrace\,t_\text{gen} / \Delta t\,\rbrace}
#'       Mean generation interval of the disease of interest
#'       in units \mjseqn{\Delta t}. FC, S, and SI methods only.
#'     }
#'     \item{`tlat`}{\mjseqn{\lbrace\,t_\text{lat} / \Delta t\,\rbrace}
#'       Mean latent period of the disease of interest
#'       in units \mjseqn{\Delta t}. SEI method only.
#'     }
#'     \item{`tinf`}{\mjseqn{\lbrace\,t_\text{inf} / \Delta t\,\rbrace}
#'       Mean infectious period of the disease of interest
#'       in units \mjseqn{\Delta t}. SEI method only.
#'     }
#'     \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'       Number of susceptible individuals at time \mjseqn{t = t_0}.
#'     }
#'     \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'       Number of exposed (infected but not infectious) individuals
#'       at time \mjseqn{t = t_0}. SEI method only.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Number of infectious individuals at time \mjseqn{t = t_0}.
#'       SI and SEI methods only.
#'     }
#'   }
#'
#' @param method A character scalar, one of `"fc"`, `"s"`, `"si"`, and `"sei"`,
#'   indicating an estimation algorithm, one of [estimate_beta_fc()],
#'   [estimate_beta_s()], [estimate_beta_si()], and [estimate_beta_sei()].
#'
#' @return
#' A fastbeta object. A list with elements:
#'
#' \describe{
#'   \item{`out`}{A data frame containing the estimation output.
#'     The result of applying [`estimate_beta_method()`][estimate-beta]
#'     to `data` and `par_list` *after* replacing missing values in
#'     `data[c("Z", "B", "mu")]` and zeros in `data$Z` with [impute_na()]
#'     (see Details).
#'   }
#'   \item{`method`}{A character scalar. The value of `method`
#'     in the function call.
#'   }
#'   \item{`warn`}{A logical scalar indicating whether a warning
#'     was issued about negative elements in `out$S` (see Details).
#'     Equal to `any(out$S < 0, na.rm = TRUE)`.
#'   }
#'   \item{`call`}{The function call. The fastbeta object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call.
#'     The fastbeta object is reproducible with `do.call(fastbeta, arg_list)`.
#'   }
#' }
#'
#' @examples
#' ## Simulate time series data using an SIR model
#' ## with seasonally forced transmission rate
#' pl <- make_par_list(model = "sir")
#' df <- make_data(pl, with_ds = TRUE, model = "sir")
#'
#' ## Estimate the seasonally forced transmission rate
#' ## using the SI method
#' fastbeta_out <- fastbeta(df, pl, method = "si")
#' plot(fastbeta_out)
#'
#' ## Simulate time series data using an SEIR model
#' ## with seasonally forced transmission rate
#' pl <- make_par_list(model = "seir")
#' df <- make_data(pl, with_ds = TRUE, model = "seir")
#'
#' ## Estimate the seasonally forced transmission rate
#' ## using the SEI method
#' fastbeta_out <- fastbeta(df, pl, method = "sei")
#' plot(fastbeta_out)
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' @seealso [methods for class "fastbeta"][fastbeta-methods],
#'   [algorithms called under the hood][estimate-beta],
#'   [try_loess()] for smoothing noisy transmission rate estimates,
#'   [bootbeta()] for constructing bootstrap confidence intervals
#'   on transmission rate estimates
#' @export
fastbeta <- function(data, par_list, method) {
  if (missing(method)) {
    stop("Missing argument `method`.")
  } else if (!(is.character(method) && length(method) == 1 &&
                 method %in% c("fc", "s", "si", "sei"))) {
    stop("`method` must be one of \"fc\", \"s\", \"si\", and \"sei\".")
  }
  if (missing(data)) {
    stop("Missing argument `data.`")
  } else if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  } else if ((method == "fc" &&
                !all(c("t", "Z", "B") %in% names(data))) ||
             (method %in% c("s", "si", "sei") &&
                !all(c("t", "Z", "B") %in% names(data)))) {
    stop("`data` is missing necessary columns.")
  }
  if (missing(par_list)) {
    stop("Missing argument `par_list`.")
  } else if (!is.list(par_list)) {
    stop("`par_list` must be a list.")
  } else if ((method %in% c("fc", "s") &&
                !all(c("S0", "tgen") %in% names(par_list))) ||
             (method == "si" &&
                !all(c("S0", "I0", "tgen") %in% names(par_list))) ||
             (method == "sei" &&
                !all(c("S0", "E0", "I0", "tlat", "tinf") %in% names(par_list)))) {
    stop("`par_list` is missing necessary elements.")
  }

  ## Save arguments in a list
  arg_list <- as.list(environment())

  ## Missing values are not tolerated. Zeros in incidence
  ## *are* tolerated but can have undesired effects.
  data[c("Z", "B", "mu")] <- mapply(impute_na,
    x = data[c("Z", "B", "mu")],
    zero_as_na = c(TRUE, FALSE, FALSE)
  )

  ## Apply the desired method
  f <- get(paste0("estimate_beta_", method))
  out <- f(data, par_list)

  ## Negative susceptibles indicates that incidence was
  ## overestimated or births were underestimated or both
  warn <- FALSE
  if (any(out$S < 0, na.rm = TRUE)) {
    warning(
      "`out$S` contains negative elements. ",
      "To prevent negative susceptibles, retry with:",
      "\n* scaled down incidence (`data$Z`), and/or",
      "\n* scaled up births (`data$B`)."
    )
    warn <- TRUE
  }

  out <- list(
    out      = out,
    method   = method,
    warn     = warn,
    call     = match.call(),
    arg_list = arg_list
  )
  structure(out, class = c("fastbeta", "data.frame"))
}
