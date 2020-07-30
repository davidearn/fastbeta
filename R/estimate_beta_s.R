#' \loadmathjax
#' Estimate time-varying transmission rates (S method)
#'
#' @description
#' `estimate_beta_s()` implements the S method (see Algorithm) for
#' estimating time-varying transmission rates \mjseqn{\beta(t)} from
#' time series of incidence, births, and natural mortality observed
#' at equally spaced time points \mjseqn{t_k = t_0 + k \Delta t}.
#' The S method is deprecated and outperformed by the more robust SI
#' method. See [estimate_beta_si()].
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
#' reproducible with `eval(call)` or `do.call(estimate_beta_s, arg_list)`.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#' The susceptible population size is estimated recursively starting from
#' the supplied initial value \mjseqn{S_0} (`par_list$S0`):
#'
#' \mjsdeqn{S_{k+1} = S_k + B_{k+1} - Z_{k+1} - \mu_k S_k \Delta t\,.}
#'
#' The infected population size is approximated by a scaling of incidence:
#' 
#' \mjsdeqn{I_k = \frac{Z_{k-g+1}}{(\gamma + \mu_k) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,,}
#' 
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_k = \frac{Z_{k+1}}{S_k I_k \Delta t}\,.}
#' 
#' ## 2. Missing data
#' 
#' Missing values in `df$Z`, `df$B`, `df$mu` are not tolerated and must
#' be imputed. Columns `S` and `beta` in the output will be filled with
#' `NA` after the index of the first missing value.
#'
#' Zeros in `df$Z` can cause divide-by-zero errors. In this case,
#' column `beta` in the output may contain `NaN` and `Inf` where
#' an estimate would otherwise be expected. Zeros in `df$Z` can be
#' "imputed" (replaced with positive values in a sensible way) to
#' avoid this issue.
#'
#' Imputation is not carried out by `estimate_beta_s()`.
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
#' df_s <- estimate_beta_s(df, par_list)
#' head(df_s)
#'
#' # Inspect
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 58) * 1e03)
#' lines(S ~ I(t - t[1]), df_s, col = "red")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.95, 1.25) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_s, col = "red")
#' 
#' @references
#' \insertRef{deJo+20}{fastbeta}
#'
#' @export
estimate_beta_s <- function(df       = data.frame(),
                            par_list = list()) {
  ## 1. Set-up -----------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "tgen")], envir = environment())

  # Preallocate memory for output
  df <- df[c("t", "Z", "B", "mu")]
  df[c("S", "I", "beta")] <- NA


  ## 2. Estimate susceptibles, infecteds, transmission rate --------------

  tgenr <- round(tgen)
  gamma <- 1 / tgen

  df[c("S", "I", "beta")] <- with(df,
    {
      S[1] <- S0
      for (i in 2:nrow(df)) {
        S[i] <- (1 - mu[i-1]) * S[i-1] + B[i] - Z[i]
      }
      if (tgenr > 0) {
        I <- c(rep(NA, tgenr - 1), Z[1:(nrow(df)-tgenr+1)]) /
          ((gamma + mu) * 1)
        beta <- (gamma + mu) * c(Z[-1], NA) /
          (S * c(rep(NA, tgenr - 1), Z[1:(nrow(df)-tgenr+1)]))
      } else {
        I <- c(Z[-1], NA) / (gamma + mu)
        beta <- (gamma + mu) * c(Z[-1], NA) /
          (S * c(Z[-1], NA))
      }
      list(S, I, beta)
    }
  )

  attr(df, "call") <- match.call()
  attr(df, "arg_list") <- arg_list
  df
}
