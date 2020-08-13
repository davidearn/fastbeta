#' \loadmathjax
#' Estimate time-varying transmission rates (SI method)
#'
#' @description
#' `estimate_beta_si()` implements the SI method (see Algorithm) for
#' estimating time-varying transmission rates \mjseqn{\beta(t)} from
#' time series of incidence, births, and natural mortality observed
#' at equally spaced time points \mjseqn{t_k = t_0 + k \Delta t}.
#'
#' @param df A data frame with numeric columns:
#'
#'   \describe{
#'     \item{`t`}{Time. `t[i]` is equal to
#'       \mjseqn{t_{i-1} = t_0 + (i-1) \Delta t} in any convenient units,
#'       where \mjseqn{t_0} is the initial observation time
#'       and \mjseqn{\Delta t} is the observation interval.
#'     }
#'     \item{`Z`}{Incidence. `Z[i]` is the number of infections between
#'       times `t[i-1]` and `t[i]`.
#'     }
#'     \item{`B`}{Births. `B[i]` is the number of births between times
#'       `t[i-1]` and `t[i]`.
#'     }
#'     \item{`mu`}{Natural mortality rate. `mu[i]` is the rate at time
#'       `t[i]` expressed per unit \mjseqn{\Delta t} and per capita.
#'     }
#'   }
#'
#' @param par_list A list with elements:
#'
#'   \describe{
#'     \item{`S0`}{Number of susceptibles at the initial observation time
#'       \mjseqn{t_0}.}
#'     \item{`I0`}{Number of infecteds at the initial observation time
#'       \mjseqn{t_0}.}
#'     \item{`tgen`}{Mean generation interval of the disease of interest
#'       in units of the observation interval \mjseqn{\Delta t}.
#'     }
#'   }
#'
#' @return
#' A data frame with numeric columns:
#'
#' \describe{
#'   \item{`t`}{Time. Identical to `df$t`.}
#'   \item{`Z`}{Incidence. Identical to `df$Z`.}
#'   \item{`B`}{Births. Identical to `df$B`.}
#'   \item{`mu`}{Natural mortality rate. Identical to `df$mu`.}
#'   \item{`S`}{Number of susceptibles. `S[i]` is the estimated
#'     number of susceptibles at time `t[i]`.
#'   }
#'   \item{`I`}{Number of infecteds. `I[i]` is the estimated
#'     number of infecteds at time `t[i]`.
#'   }
#'   \item{`beta`}{Transmission rate. `beta[i]` is the estimated
#'     transmission rate at time `t[i]`, expressed per unit \mjseqn{\Delta t}
#'     per susceptible per infected.
#'   }
#' }
#'
#' The data frame has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(estimate_beta_si, arg_list)`.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#' The susceptible and infected population sizes are estimated recursively
#' starting from the supplied initial values \mjseqn{S_0} (`par_list$S0`)
#' and \mjseqn{I_0} (`par_list$I0`):
#'
#' \mjsdeqn{\begin{align*} S_{k+1} &= \frac{\big\lbrack 1 - \frac{1}{2} \mu_k \Delta t \big\rbrack S_k + B_{k+1} - Z_{k+1}}{1 + \frac{1}{2} \mu_{k+1} \Delta t}\,, \cr I_{k+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\gamma + \mu_k) \Delta t \big\rbrack I_k + Z_{k+1}}{1 + \frac{1}{2} (\gamma + \mu_{k+1}) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,. \end{align*}}
#'
#' The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_k = \frac{Z_k + Z_{k+1}}{2 S_k I_k \Delta t}\,.}
#' 
#' ## 2. Missing data
#' 
#' Missing values in `df$Z`, `df$B`, `df$mu` are not tolerated and
#' must be imputed. Columns `S`, `I`, `beta` in the output will be
#' filled with `NA` after the index of the first missing value.
#'
#' Strings of zeros in `df$Z` can lead to spurious zeros and large
#' numeric elements in column `beta` of the output. Zeros in `df$Z`
#' can be "imputed" (replaced with positive values in a sensible way)
#' to avoid this issue.
#'
#' Imputation is not carried out by `estimate_beta_si()`.
#'
#' @examples
#' # Simulate time series of incidence and births using
#' # an SIR model with seasonally forced transmission rate
#' par_list <- make_par_list(dt_weeks = 1, epsilon = 0.5)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7, # 20 years is ~1042 weeks
#'   with_dem_stoch = TRUE
#' )
#' head(df)
#'
#' # Estimate susceptibles, infecteds, and
#' # the seasonally forced transmission rate
#' df_si <- estimate_beta_si(df, par_list)
#' head(df_si)
#'
#' # Inspect
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 58) * 1e03)
#' lines(S ~ I(t - t[1]), df_si, col = "red")
#' plot(I ~ I(t - t[1]), df, type = "l", ylim = c(3, 3) * 1e03)
#' lines(I ~ I(t - t[1]), df_si, col = "red")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.95, 1.25) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_si, col = "red")
#'
#' @references
#' \insertRef{deJo+20}{fastbeta}
#'
#' @export
estimate_beta_si <- function(df       = data.frame(),
                             par_list = list()) {
  ## 1. Set-up -----------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "I0", "tgen")], envir = environment())

  # Preallocate memory for output
  df <- df[c("t", "Z", "B", "mu")]
  df[c("S", "I", "beta")] <- NA


  ## 2. Estimate susceptibles, infecteds, transmission rate --------------

  gamma <- 1 / tgen

  df[c("S", "I", "beta")] <- with(df,
    {
      S[1] <- S0
      I[1] <- I0
      for (i in 2:nrow(df)) {
        S[i] <- ((1 - 0.5 * mu[i-1] * 1) * S[i-1] + B[i] - Z[i]) /
          (1 + 0.5 * mu[i] * 1)
        I[i] <- ((1 - 0.5 * (gamma + mu[i-1]) * 1) * I[i-1] + Z[i]) /
          (1 + 0.5 * (gamma + mu[i]) * 1)
      }
      beta <- (Z + c(Z[-1], NA)) / (2 * S * I * 1)
      list(S, I, beta)
    }
  )

  attr(df, "call") <- match.call()
  attr(df, "arg_list") <- arg_list
  df
}
