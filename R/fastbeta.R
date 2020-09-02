#' \loadmathjax
#' Estimate time-varying transmission rates
#'
#' @description
#' A wrapper for the [`estimate_beta_method()`][estimate-beta] functions,
#' which implement the FC, S, SI, and SEI methods for estimating time-varying
#' transmission rates \mjseqn{\beta(t)} from equally spaced time series of
#' incidence, births, and natural mortality.
#'
#' @details
#' The help page for [`estimate_beta_method()`][estimate-beta]
#' documents the four estimation algorithms, their output,
#' various issues, and possible fixes. To summarize the main
#' issues, the algorithms struggle with:
#'
#' 1. Missing data and zeros in incidence.
#' 2. Systematic overestimation of incidence by `data$Z`
#'    and systematic underestimation of births by `data$B`,
#'    causing the estimated susceptible population size
#'    to become negative.
#' 3. Process and observation error in the data-generating
#'    process, causing spurious noise in the estimated
#'    transmission rate.
#'
#' `fastbeta()` handles these issues as follows:
#'
#' 1. Prior to calling [`estimate_beta_method()`][estimate-beta],
#'    `fastbeta()` uses [impute_na()] to replace missing values in
#'    `data[c("Z", "B", "mu")]` and zeros in `data$Z` with positive
#'    values. [impute_na()], which performs linear interpolation,
#'    fails in certain cases that are documented in its help page.
#'    If a more sophisticated technique is desired, then it must be
#'    applied prior to the call to `fastbeta()`.
#' 2. `fastbeta()` warns if the estimated susceptible population size
#'    is at any time negative and suggests restarting with scaled down
#'    `data$Z` and/or scaled up `data$B`. (Less crude approaches to
#'    correcting incidence certainly exist. For an approach based on
#'    deconvolution, see [deconvol()].) A rule-of-thumb is to iterate
#'    until the start of the susceptible time series no longer displays
#'    pronounced transient dynamics, while still considering whether
#'    the correction applied is reasonable.
#' 3. `fastbeta()` *does not* prevent or smooth out noise in the estimated
#'    transmission rate. Determining what constitutes a "good" degree of
#'    smoothing is not a trivial task and typically involves a graphical
#'    inspection of results for different values of smoothing parameters.
#'    [try_loess()] attempts to facilitate this process for users fitting
#'    smooth loess curves to noisy time series.
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
#' # Simulate time series data using an SIR model
#' # with seasonally forced transmission rate
#' pl <- make_par_list(model = "sir")
#' df <- make_data(pl, with_ds = TRUE, model = "sir")
#'
#' # Estimate the seasonally forced transmission rate
#' # using the SI method
#' fastbeta_out <- fastbeta(df, pl, method = "si")
#' plot(fastbeta_out)
#'
#' # Simulate time series data using an SEIR model
#' # with seasonally forced transmission rate
#' pl <- make_par_list(model = "seir")
#' df <- make_data(pl, with_ds = TRUE, model = "seir")
#'
#' # Estimate the seasonally forced transmission rate
#' # using the SEI method
#' fastbeta_out <- fastbeta(df, pl, method = "sei")
#' plot(fastbeta_out)
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' @seealso [methods for class "fastbeta"][fastbeta-methods],
#'   [algorithms used under the hood][estimate-beta],
#'   [try_loess()] for smoothing noisy transmission rate estimates,
#'   [bsbeta()] for constructing bootstrap confidence intervals on
#'   transmission rate estimates
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

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Missing values are not tolerated. Zeros in incidence
  # *are* tolerated but can have undesired effects.
  data[c("Z", "B", "mu")] <- mapply(impute_na,
    x = data[c("Z", "B", "mu")],
    zero_as_na = c(TRUE, FALSE, FALSE)
  )

  # Apply the desired method
  f <- get(paste0("estimate_beta_", method))
  out <- f(data, par_list)

  # Negative susceptibles indicates that incidence was
  # overestimated or births were underestimated or both
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
