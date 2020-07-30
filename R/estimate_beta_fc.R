#' \loadmathjax
#' Estimate time-varying transmission rates (FC method)
#'
#' @description
#' `estimate_beta_fc()` implements the FC method (see Algorithm) for
#' estimating time-varying transmission rates \mjseqn{\beta(t)} from
#' time series of incidence and births observed at equally spaced
#' time points \mjseqn{t_k = t_0 + k \Delta t}. The FC method is
#' deprecated and outperformed by the more robust SI method. See
#' [estimate_beta_si()].
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
#'   \item{`Z_agg`}{Aggregated incidence. `Z_agg[i]` is the number
#'     of infections between times `t[i-round(tgen)+1]` and `t[i]`,
#'     for `i` in `seq(1 + round(tgen), nrow(df), by = round(tgen))`.
#'     `Z_agg[i]` is `NA` otherwise.
#'   }
#'   \item{`B`}{Births. Identical to `df$B`.}
#'   \item{`B_agg`}{Aggregated births. `B_agg[i]` is the number
#'     of births between times `t[i-round(tgen)+1]` and `t[i]`,
#'     for `i` in `seq(1 + round(tgen), nrow(df), by = round(tgen))`.
#'     `B_agg[i]` is `NA` otherwise.
#'   }
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
#' reproducible with `eval(call)` or `do.call(estimate_beta_fc, arg_list)`.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#' The supplied incidence time series \mjseqn{Z_k} (`df$Z`) is aggregated
#' over generation intervals:
#'
#' \mjsdeqn{Z_k^\text{agg} = \sum_{i=k-g+1}^k Z_i\,,\quad k=jg\,,\quad j = 1,2,\ldots\,,}
#'
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The supplied births time series \mjseqn{B_k} (`df$B`)
#' is aggregated similarly. The susceptible population size is estimated
#' recursively starting from the supplied initial value \mjseqn{S_0}
#' (`par_list$S0`):
#'
#' \mjsdeqn{S_{k+g} = S_k + B_{k+g}^\text{agg} - Z_{k+g}^\text{agg}\,,\quad k=jg\,,\quad j = 0,1,\ldots\,,}
#'
#' where natural mortality of susceptibles is assumed to be negligible.
#' The infected population size is approximated by aggregated incidence:
#'
#' \mjsdeqn{I_k = Z_k^\text{agg}\,,\quad k=jg\,,\quad j = 1,2,\ldots\,,}
#' 
#' where natural mortality of infecteds is assumed to be negligible.
#' The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_k = \frac{Z_{k+g}}{S_k I_k g \Delta t}\,,\quad k=jg\,,\quad j = 1,2,\ldots\,.}
#' 
#' ## 2. Missing data
#' 
#' Missing values in `df$Z` and `df$B` are not tolerated and must be
#' imputed. Columns `S` and `beta` in the output will be filled with
#' `NA` after the index of the first missing value.
#'
#' Zeros in `df$Z` can cause divide-by-zero errors. In this case,
#' column `beta` in the output may contain `NaN` and `Inf` where
#' an estimate would otherwise be expected. Zeros in `df$Z` can be
#' "imputed" (replaced with positive values in a sensible way) to
#' avoid this issue.
#'
#' Imputation is not carried out by `estimate_beta_fc()`.
#'
#' ## 3. Effective observation interval
#' 
#' In the absence of errors, every `round(par_list$tgen)`th row in
#' the output will be complete. The remaining rows will contain `NA`
#' in columns `S`, `I`, `beta`, `Z_agg`, and `B_agg`. This is not a
#' bug: the (rounded) generation interval becomes the *effective*
#' observation interval when the incidence and birth time series
#' are aggregated (see Algorithm).
#'
#' @examples
#' # Simulate time series of incidence and births using
#' # an SIR model with seasonally forced transmission rate
#' par_list <- make_par_list(dt_weeks = 1)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7, # 20 years is ~1042 weeks
#'   with_dem_stoch = FALSE
#' )
#' head(df)
#'
#' # Estimate susceptibles, infecteds, and
#' # the seasonally forced transmission rate
#' df_fc <- estimate_beta_fc(df, par_list)
#' head(df_fc)
#'
#' # Estimation of susceptibles and transmission rate fails
#' # because the FC method ignores susceptible mortality
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 83) * 1e03)
#' lines(S ~ I(t - t[1]), df_fc[!is.na(df_fc$S), ], col = "red")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.5, 1.2) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_fc[!is.na(df_fc$beta), ], col = "red")
#' 
#' @references
#' \insertRef{deJo+20}{fastbeta}
#'
#' @export
estimate_beta_fc <- function(df       = data.frame(),
                             par_list = list()) {
  ## 1. Set-up -----------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "tgen")], envir = environment())

  # Preallocate memory for output
  df <- data.frame(
    t     = df$t,
    Z     = df$Z,
    Z_agg = NA,
    B     = df$B,
    B_agg = NA,
    S     = NA,
    I     = NA,
    beta  = NA
  )


  ## 2. Aggregate incidence, births --------------------------------------

  # Aggregates are recorded after each generation interval,
  # starting at the first time point
  tgenr <- round(tgen)
  ind_with_agg <- seq(1, nrow(df), by = tgenr)

  # Aggregate at the first time point cannot be computed,
  # because there are no prior observations to aggregate
  for (i in ind_with_agg[-1]) {
    ind_into_agg <- seq(i - tgenr + 1, i, by = 1)
    df$B_agg[i] <- sum(df$B[ind_into_agg])
    df$Z_agg[i] <- sum(df$Z[ind_into_agg])
  }


  ## 3. Estimate susceptibles, infecteds, transmission rate --------------

  df[c("S", "I", "beta")] <- with(df,
    {
      # 
      S[1] <- S0
      for (i in ind_with_agg[-1]) {
        S[i] <- S[i-tgenr] + B_agg[i] - Z_agg[i]
      }
      I <- Z_agg
      beta[ind_with_agg] <- c(Z_agg[ind_with_agg[-1]], NA) /
        (S[ind_with_agg] * Z_agg[ind_with_agg] * tgenr)
      list(S, I, beta)
    }
  )

  attr(df, "call") <- match.call()
  attr(df, "arg_list") <- arg_list
  df
}
