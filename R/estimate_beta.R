#' \loadmathjax
#' Estimate time-varying transmission rates
#'
#' @description
#' Functions implementing the FC, S, and SI methods (see Algorithm)
#' for estimating time-varying transmission rates \mjseqn{\beta(t)}
#' from time series of incidence, births, and natural mortality with
#' observations at equally spaced time points
#' \mjseqn{t_i = t_0 + i \Delta t}.
#' The FC and S methods are deprecated and outperformed by the
#' more robust SI method. `estimate_beta_si()` need not be called
#' directly: it is wrapped in the more useful constructor function
#' [fastbeta()], which does checks on input and output and returns
#' a fastbeta object with associated methods.
#' 
#' @details
#' # Details
#'
#' ## 1. Algorithm
#'
#' ### FC method
#' The supplied incidence time series \mjseqn{Z_i} (`df$Z`) is aggregated
#' over generation intervals:
#'
#' \mjsdeqn{Z_i^\text{agg} = \sum_{k=i-g+1}^k Z_k\,,\quad i=jg\,,\quad j = 1,2,\ldots\,,}
#'
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The supplied births time series \mjseqn{B_i} (`df$B`)
#' is aggregated similarly. The susceptible population size is estimated
#' recursively starting from the supplied initial value \mjseqn{S_0}
#' (`par_list$S0`):
#'
#' \mjsdeqn{S_{i+g} = S_i + B_{i+g}^\text{agg} - Z_{i+g}^\text{agg}\,,\quad i=jg\,,\quad j = 0,1,\ldots\,,}
#'
#' where natural mortality of susceptibles is assumed to be negligible.
#' The infected population size is approximated by aggregated incidence:
#'
#' \mjsdeqn{I_i = Z_i^\text{agg}\,,\quad i=jg\,,\quad j = 1,2,\ldots\,,}
#' 
#' where natural mortality of infecteds is assumed to be negligible.
#' The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_{i+g}}{S_i I_i g \Delta t}\,,\quad i=jg\,,\quad j = 1,2,\ldots\,.}
#'
#' Note that this algorithm was first published by
#' \insertCite{FineClar82;textual}{fastbeta}.
#'
#' ### S method
#' The susceptible population size is estimated recursively starting from
#' the supplied initial value \mjseqn{S_0} (`par_list$S0`):
#'
#' \mjsdeqn{S_{i+1} = S_i + B_{i+1} - Z_{i+1} - \mu_i S_i \Delta t\,.}
#'
#' The infected population size is approximated by a scaling of incidence:
#' 
#' \mjsdeqn{I_i = \frac{Z_{i-g+1}}{(\gamma + \mu_i) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,,}
#' 
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_{i+1}}{S_i I_i \Delta t}\,.}
#'
#' ### SI method
#' The susceptible and infected population sizes are estimated recursively
#' starting from the supplied initial values \mjseqn{S_0} (`par_list$S0`)
#' and \mjseqn{I_0} (`par_list$I0`):
#'
#' \mjsdeqn{\begin{align*} S_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} \mu_i \Delta t \big\rbrack S_i + B_{i+1} - Z_{i+1}}{1 + \frac{1}{2} \mu_{i+1} \Delta t}\,, \cr I_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\gamma + \mu_i) \Delta t \big\rbrack I_i + Z_{i+1}}{1 + \frac{1}{2} (\gamma + \mu_{i+1}) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,. \end{align*}}
#'
#' The transmission rate is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_i + Z_{i+1}}{2 S_i I_i \Delta t}\,.}
#' 
#' ## 2. Missing data
#' 
#' Missing values in `df$Z` are not tolerated and must be imputed.
#' Columns `S` and `beta` in the output will be filled with `NA`
#' after the index of the first missing value.
#'
#' In the FC and S methods, zeros in `df$Z` can cause divide-by-zero
#' errors. In this case, column `beta` in the output may contain `NaN`
#' and `Inf` where an estimate would otherwise be expected. In the SI
#' method, strings of zeros in `df$Z` can lead to spurious zeros and
#' large numeric elements in column `beta` of the output. Zeros in `df$Z`
#' can be "imputed" (replaced with positive values in a sensible way)
#' to avoid this issue.
#'
#' Imputation is not carried out by `estimate_beta_fc()` and friends,
#' but *is* carried out by the more useful constructor function
#' [fastbeta()], a wrapper for `estimate_beta_si()`.
#'
#' ## 3. Effective observation interval (FC method only)
#' 
#' In the absence of errors, every `round(par_list$tgen)`th row in
#' the FC method output will be complete. The remaining rows will
#' contain `NA` in columns `S`, `I`, `beta`, `Z_agg`, and `B_agg`.
#' This is not a mistake: the rounded generation interval becomes
#' the *effective* observation interval when the incidence and births
#' time series are aggregated (see Algorithm).
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
#'       S and SI methods only.
#'     }
#'   }
#' 
#' @param par_list A list with elements:
#'
#'   \describe{
#'     \item{`tgen`}{Mean generation interval of the disease of interest
#'       in units of the observation interval \mjseqn{\Delta t}.
#'     }
#'     \item{`S0`}{Number of susceptibles at the initial observation time
#'       \mjseqn{t_0}.
#'     }
#'     \item{`I0`}{Number of infecteds at the initial observation time
#'       \mjseqn{t_0}. SI method only.
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
#'     `Z_agg[i]` is `NA` otherwise. FC method only.
#'   }
#'   \item{`B`}{Births. Identical to `df$B`.}
#'   \item{`B_agg`}{Aggregated births. `B_agg[i]` is the number
#'     of births between times `t[i-round(tgen)+1]` and `t[i]`,
#'     for `i` in `seq(1 + round(tgen), nrow(df), by = round(tgen))`.
#'     `B_agg[i]` is `NA` otherwise. FC method only.
#'   }
#'   \item{`mu`}{Natural mortality rate. Identical to `df$mu`.
#'     S and SI methods only.
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
#' reproducible with `eval(call)` or `do.call(estimate_beta_x, arg_list)`.
#' 
#' @examples
#' # Simulate time series data using an SIR model
#' # with seasonally forced transmission rate
#' par_list <- make_par_list(dt_weeks = 1)
#' df <- make_data(par_list = par_list, with_ds = FALSE)
#' head(df)
#'
#' # Estimate susceptibles, infecteds, and
#' # the seasonally forced transmission rate
#' df_fc <- estimate_beta_fc(df, par_list)
#' head(df_fc)
#' df_s <- estimate_beta_s(df, par_list)
#' head(df_s)
#' df_si <- estimate_beta_si(df, par_list)
#' head(df_si)
#'
#' # FC method fails rapidly as a result
#' # of ignoring susceptible mortality
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 83) * 1e03)
#' lines(S ~ I(t - t[1]), df_fc[!is.na(df_fc$S), ], col = "blue")
#' plot(I ~ I(t - t[1]), df, type = "l", ylim = c(0, 3) * 1e03)
#' lines(I ~ I(t - t[1]), df_fc[!is.na(df_fc$I), ], col = "blue")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.5, 1.2) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_fc[!is.na(df_fc$beta), ], col = "blue")
#'
#' # S method does well
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 58) * 1e03)
#' lines(S ~ I(t - t[1]), df_s, col = "blue")
#' plot(I ~ I(t - t[1]), df, type = "l", ylim = c(0, 3) * 1e03)
#' lines(I ~ I(t - t[1]), df_s, col = "red")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.95, 1.25) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_s, col = "blue")
#'
#' # SI method does well
#' plot(S ~ I(t - t[1]), df, type = "l", ylim = c(43, 58) * 1e03)
#' lines(S ~ I(t - t[1]), df_si, col = "blue")
#' plot(I ~ I(t - t[1]), df, type = "l", ylim = c(3, 3) * 1e03)
#' lines(I ~ I(t - t[1]), df_si, col = "blue")
#' plot(beta ~ I(t - t[1]), df, type = "l", ylim = c(0.95, 1.25) * 1e-05)
#' lines(beta ~ I(t - t[1]), df_si, col = "blue")
#'
#' # However, restarting with `with_ds = TRUE`
#' # reveals that the S method is vulnerable
#' # to propagation of noise from the data to
#' # the transmission rate estimate. Hence the
#' # SI method should be preferred in practice.
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#' \insertRef{FineClar82}{fastbeta}
#'
#' @seealso [fastbeta()]
#' @name estimate_beta
NULL

#' @rdname estimate_beta
#' @export
estimate_beta_fc <- function(df, par_list) {
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

  structure(df, call = match.call(), arg_list = arg_list)
}

#' @rdname estimate_beta
#' @export
estimate_beta_s <- function(df, par_list) {
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

  structure(df, call = match.call(), arg_list = arg_list)
}

#' @rdname estimate_beta
#' @export
estimate_beta_si <- function(df, par_list) {
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

  structure(df, call = match.call(), arg_list = arg_list)
}
