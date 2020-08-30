#' \loadmathjax
#' Estimate time-varying transmission rates
#'
#' @description
#' Functions implementing the FC, S, SI, and SEI methods (see Algorithm)
#' for estimating time-varying transmission rates \mjseqn{\beta(t)} from
#' time series of incidence, births, and natural mortality with observations
#' at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}. The FC
#' and S methods are deprecated and outperformed by the more robust SI
#' and SEI methods. These four functions should not be called directly,
#' as they are all accessible through the `method` argument of the more
#' useful constructor function [fastbeta()], which does checks on input
#' and output and returns a fastbeta object with associated methods.
#'
#' @details
#' # Details
#'
#' ## 1. Algorithm
#'
#' ### FC method
#' The supplied incidence time series \mjseqn{Z_i} (`data$Z`) is aggregated
#' over generation intervals:
#'
#' \mjsdeqn{Z_i^\text{agg} = \sum_{k=i-g+1}^k Z_k\,,\quad i=jg\,,\quad j = 1,2,\ldots\,,}
#'
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The supplied births time series \mjseqn{B_i} (`data$B`)
#' is aggregated similarly. The susceptible population size is estimated
#' recursively starting from the supplied initial value \mjseqn{S_0}
#' (`par_list$S0`):
#'
#' \mjsdeqn{S_{i+g} = S_i + B_{i+g}^\text{agg} - Z_{i+g}^\text{agg}\,,\quad i=jg\,,\quad j = 0,1,\ldots\,,}
#'
#' where natural mortality of susceptibles is assumed to be negligible.
#' The infectious population size is approximated by aggregated incidence:
#'
#' \mjsdeqn{I_i = Z_i^\text{agg}\,,\quad i=jg\,,\quad j = 1,2,\ldots\,,}
#'
#' where natural mortality of infectious individuals is assumed to
#' be negligible. The transmission rate per susceptible individual
#' per infectious individual is then estimated as
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
#' The infectious population size is approximated by a scaling of incidence:
#'
#' \mjsdeqn{I_i = \frac{Z_{i-g+1}}{(\gamma + \mu_i) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,,}
#'
#' where \mjseqn{g = \mathrm{nint}(t_\text{gen} / \Delta t)} is the mean
#' generation interval \mjseqn{t_\text{gen}} in units of the observation
#' interval \mjseqn{\Delta t} (`par_list$tgen`), rounded to the nearest
#' integer. The transmission rate per susceptible individual
#' per infectious individual is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_{i+1}}{S_i I_i \Delta t}\,.}
#'
#' ### SI method
#' The susceptible and infectious population sizes are estimated recursively
#' starting from the supplied initial values \mjseqn{S_0} (`par_list$S0`)
#' and \mjseqn{I_0} (`par_list$I0`):
#'
#' \mjsdeqn{\begin{align*} S_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} \mu_i \Delta t \big\rbrack S_i + B_{i+1} - Z_{i+1}}{1 + \frac{1}{2} \mu_{i+1} \Delta t}\,, \cr I_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\gamma + \mu_i) \Delta t \big\rbrack I_i + Z_{i+1}}{1 + \frac{1}{2} (\gamma + \mu_{i+1}) \Delta t}\,,\quad \gamma = 1 / t_\text{gen}\,. \end{align*}}
#'
#' The transmission rate per susceptible individual
#' per infectious individual is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_i + Z_{i+1}}{2 S_i I_i \Delta t}\,.}
#'
#' ### SEI method
#' The susceptible, exposed (infected but not infectious), and infectious
#' population sizes are estimated recursively starting from the supplied
#' initial values \mjseqn{S_0} (`par_list$S0`), \mjseqn{E_0} (`par_list$E0`),
#' and \mjseqn{I_0} (`par_list$I0`):
#'
#' \mjsdeqn{\begin{align*} S_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} \mu_i \Delta t \big\rbrack S_i + B_{i+1} - Z_{i+1}}{1 + \frac{1}{2} \mu_{i+1} \Delta t}\,, \cr E_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\sigma + \mu_i) \Delta t \big\rbrack E_i + Z_i}{1 + \frac{1}{2} (\sigma + \mu_{i+1}) \Delta t}\,,\quad \sigma = 1 / t_\text{lat}\,, \cr I_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\gamma + \mu_i) \Delta t \big\rbrack I_i + \frac{1}{2} \sigma (E_i + E_{i+1}) \Delta t}{1 + \frac{1}{2} (\gamma + \mu_{i+1}) \Delta t}\,,\quad \sigma = 1 / t_\text{lat}\,,\quad \gamma = 1 / t_\text{inf}\,. \end{align*}}
#'
#' The transmission rate per susceptible individual
#' per infectious individual is then estimated as
#'
#' \mjsdeqn{\beta_i = \frac{Z_i + Z_{i+1}}{2 S_i I_i \Delta t}\,.}
#'
#' ## 2. Missing data and zeros in incidence
#'
#' Missing values in `data$Z` are not tolerated and must be imputed.
#' Columns `S` and `beta` in the output will be filled with `NA`
#' after the index of the first missing value. In the FC and S methods,
#' zeros in `data$Z` can cause divide-by-zero errors. In this case,
#' column `beta` in the output may contain `NaN` and `Inf` where an
#' estimate would otherwise be expected. In the SI and SEI methods,
#' strings of zeros in `data$Z` can lead to spurious zeros and large
#' numeric elements in column `beta` of the output. Zeros in `data$Z`
#' can be "imputed" (replaced with positive values in a sensible way)
#' to avoid these issues. `estimate_beta_method()` does not attempt
#' imputation. However, see [fastbeta()].
#'
#' ## 3. Negative susceptibles
#'
#' If true incidence is systematically overestimated by `data$Z`
#' or true births are systematically underestimated by `data$B`,
#' then the estimated susceptible population size can grossly
#' underestimate the true number of susceptibles and can even
#' become negative. In the latter case, column `S` in the output
#' will contain negative elements. In both cases, more reasonable
#' results are typically obtained by restarting with scaled up
#' `data$Z` and/or scaled down `data$B`. `estimate_beta_method()`
#' does not warn about negative susceptibles. However, see [fastbeta()].
#'
#' ## 4. Noise in the estimated transmission rate
#'
#' All four estimation methods propagate noise (due to process
#' and observation error) from the data to transmission rate
#' estimates. To distill temporal patterns of interest from the
#' noise, the raw transmission rate estimates can be smoothed
#' using, for example, [stats::loess()]. `estimate_beta_method()`
#' does not undertake smoothing.
#'
#' ## 5. Effective observation interval (FC method only)
#'
#' In the absence of errors, every `round(par_list$tgen)`th row in
#' the FC method output will be complete. The remaining rows will
#' contain `NA` in columns `S`, `I`, `beta`, `Z_agg`, and `B_agg`.
#' This is not a mistake: the rounded generation interval becomes
#' the *effective* observation interval when the incidence and births
#' time series are aggregated (see Algorithm).
#'
#' @inheritParams fastbeta
#'
#' @return
#' A data frame with numeric columns:
#'
#' \describe{
#'   \item{`t`}{\mjseqn{\lbrace\,t_i\,\rbrace}
#'     Time. Identical to `data$t`.
#'   }
#'   \item{`Z`}{\mjseqn{\lbrace\,Z_i\,\rbrace}
#'     Incidence. Identical to `data$Z`.
#'   }
#'   \item{`Z_agg`}{\mjseqn{\lbrace\,Z_i^\text{agg}\,\rbrace}
#'     Aggregated incidence. `Z_agg[i]` is the number of infections
#'     between times `t[i-round(tgen)+1]` and `t[i]`, for `i` in
#'     `seq(1 + round(tgen), nrow(data), by = round(tgen))`.
#'     `Z_agg[i]` is `NA` for all other `i`. FC method only.
#'   }
#'   \item{`B`}{\mjseqn{\lbrace\,B_i\,\rbrace}
#'     Births. Identical to `data$B`.
#'   }
#'   \item{`B_agg`}{\mjseqn{\lbrace\,B_i^\text{agg}\,\rbrace}
#'     Aggregated births. `B_agg[i]` is the number of births
#'     between times `t[i-round(tgen)+1]` and `t[i]`, for `i` in
#'     `seq(1 + round(tgen), nrow(data), by = round(tgen))`.
#'     `B_agg[i]` is `NA` for all other `i`. FC method only.
#'   }
#'   \item{`mu`}{\mjseqn{\lbrace\,\mu_i \Delta t\,\rbrace}
#'     Per capita natural mortality rate. Identical to `data$mu`.
#'     S, SI, and SEI methods only.
#'   }
#'   \item{`S`}{\mjseqn{\lbrace\,S_i\,\rbrace}
#'     Estimated number of susceptible individuals.
#'   }
#'   \item{`E`}{\mjseqn{\lbrace\,E_i\,\rbrace}
#'     Estimated number of exposed (infected but not infectious) individuals.
#'     SEI method only.
#'   }
#'   \item{`I`}{\mjseqn{\lbrace\,I_i\,\rbrace}
#'     Estimated number of infectious individuals.
#'   }
#'   \item{`beta`}{\mjseqn{\lbrace\,\beta_i \Delta t\,\rbrace}
#'     Estimated transmission rate, expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'   }
#' }
#'
#' @examples
#' # Simulate time series data using an SIR model
#' # to test the FC, S, and SI methods
#' pl_sir <- make_par_list(model = "sir")
#' df_sir <- make_data(pl_sir, with_ds = FALSE, model = "sir")
#' head(df_sir)
#'
#' # Simulate time series data using an SEIR model
#' # to test the SEI method
#' pl_seir <- make_par_list(model = "seir")
#' df_seir <- make_data(pl_seir, with_ds = FALSE, model = "seir")
#' head(df_seir)
#'
#' # Estimate the seasonally forced transmission rate
#' # underlying each simulation
#' f_list <- list(
#'   fc  = estimate_beta_fc,
#'   s   = estimate_beta_s,
#'   si  = estimate_beta_si,
#'   sei = estimate_beta_sei
#' )
#' eb_out <- mapply(function(f, df, pl) f(data = df, par_list = pl),
#'   f = f_list,
#'   df = list(df_sir, df_sir, df_sir, df_seir),
#'   pl = list(pl_sir, pl_sir, pl_sir, pl_seir),
#'   SIMPLIFY = FALSE
#' )
#' lapply(eb_out, head)
#'
#' op <- par(mfrow = c(3, 1), mar = c(2, 6, 1, 1), oma = c(3, 0, 2, 0))
#'
#' # FC method fails rapidly as a result of ignoring
#' # susceptible mortality
#' plot(S ~ t, df_sir, type = "l", lwd = 4, ylim = c(43, 83) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, na.omit(eb_out$fc), lwd = 2, col = "seagreen")
#' title(main = "FC method", line = 1, xpd = NA)
#' plot(I ~ t, df_sir, type = "l", lwd = 4, ylim = c(0, 3.5) * 1e03,
#'      xlab = "", ylab = "Infectious")
#' lines(I ~ t, na.omit(eb_out$fc), lwd = 2, col = "mediumvioletred")
#' plot(beta ~ t, df_sir, type = "l", lwd = 4, ylim = c(0.5, 1.3) * 1e-05,
#'      xlab = "", ylab = "Transmission rate")
#' lines(beta ~ t, na.omit(eb_out$fc), lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' # S method is visibly biased
#' plot(S ~ t, df_sir, type = "l", lwd = 4, ylim = c(43, 58) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, eb_out$s, lwd = 2, col = "seagreen")
#' title(main = "S method", line = 1, xpd = NA)
#' plot(I ~ t, df_sir, type = "l", lwd = 4, ylim = c(0, 3.5) * 1e03,
#'      xlab = "", ylab = "Infectious")
#' lines(I ~ t, eb_out$s, lwd = 2, col = "mediumvioletred")
#' plot(I(beta / pl_sir$beta_mean) ~ t, df_sir, type = "l", lwd = 4, ylim = c(0.7, 1.3),
#'      xlab = "", ylab = "Transmission rate\nrelative to mean")
#' lines(I(beta / pl_sir$beta_mean) ~ t, eb_out$s, lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' # SI method does well
#' plot(S ~ t, df_sir, type = "l", lwd = 4, ylim = c(43, 58) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, eb_out$si, lwd = 2, col = "seagreen")
#' title(main = "SI method", line = 1, xpd = NA)
#' plot(I ~ t, df_sir, type = "l", lwd = 4, ylim = c(0, 3.5) * 1e03,
#'      xlab = "", ylab = "Infectious")
#' lines(I ~ t, eb_out$si, lwd = 2, col = "mediumvioletred")
#' plot(I(beta / pl_sir$beta_mean) ~ t, df_sir, type = "l", lwd = 4, ylim = c(0.7, 1.3),
#'      xlab = "", ylab = "Transmission rate\nrelative to mean")
#' lines(I(beta / pl_sir$beta_mean) ~ t, eb_out$si, lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' # SEI method does well
#' plot(S ~ t, df_seir, type = "l", lwd = 4, ylim = c(43, 58) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, eb_out$sei, lwd = 2, col = "seagreen")
#' title(main = "SEI method", line = 1, xpd = NA)
#' plot(I ~ t, df_seir, type = "l", lwd = 4, ylim = c(0, 2) * 1e03,
#'      xlab = "", ylab = "Infectious")
#' lines(I ~ t, eb_out$sei, lwd = 2, col = "mediumvioletred")
#' plot(I(beta / pl_seir$beta_mean) ~ t, df_seir, type = "l", lwd = 4, ylim = c(0.7, 1.3),
#'      xlab = "", ylab = "Transmission rate\nrelative to mean")
#' lines(I(beta / pl_seir$beta_mean) ~ t, eb_out$sei, lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' par(op)
#'
#' # NOTE: Restarting with `with_ds = TRUE` passed to
#' # `make_data()` suggests that the SI method is the
#' # most robust to demographic stochasticity and more
#' # generally to noise in the data-generating process.
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#' \insertRef{FineClar82}{fastbeta}
#'
#' @seealso [fastbeta()]
#' @name estimate-beta
NULL

#' @rdname estimate-beta
#' @export
estimate_beta_fc <- function(data, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "tgen")], envir = environment())

  # Preallocate memory for output
  data <- data.frame(
    t     = data$t,
    Z     = data$Z,
    Z_agg = NA,
    B     = data$B,
    B_agg = NA,
    S     = NA,
    I     = NA,
    beta  = NA
  )


  ## Aggregate incidence, births -----------------------------------------

  # Aggregates are recorded after each generation interval,
  # starting at the first time point
  tgenr <- round(tgen)
  ind_with_agg <- seq(1, nrow(data), by = tgenr)

  # Aggregate at the first time point cannot be computed,
  # because there are no prior observations to aggregate
  for (i in ind_with_agg[-1]) {
    ind_into_agg <- seq(i - tgenr + 1, i, by = 1)
    data$B_agg[i] <- sum(data$B[ind_into_agg])
    data$Z_agg[i] <- sum(data$Z[ind_into_agg])
  }


  ## Estimate susceptible, infectious population sizes ...  -------------
  ## and transmission rate

  data[c("S", "I", "beta")] <- with(data,
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

  data
}

#' @rdname estimate-beta
#' @export
estimate_beta_s <- function(data, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "tgen")], envir = environment())

  # Preallocate memory for output
  data <- data[c("t", "Z", "B", "mu")]
  data[c("S", "I", "beta")] <- NA


  ## Estimate susceptible, infectious population sizes ...  -------------
  ## and transmission rate

  tgenr <- round(tgen)
  gamma <- 1 / tgen

  data[c("S", "I", "beta")] <- with(data,
    {
      S[1] <- S0
      for (i in 2:nrow(data)) {
        S[i] <- (1 - mu[i-1]) * S[i-1] + B[i] - Z[i]
      }
      if (tgenr > 0) {
        I <- c(rep(NA, tgenr - 1), Z[1:(nrow(data)-tgenr+1)]) /
          ((gamma + mu) * 1)
        beta <- (gamma + mu) * c(Z[-1], NA) /
          (S * c(rep(NA, tgenr - 1), Z[1:(nrow(data)-tgenr+1)]))
      } else {
        I <- c(Z[-1], NA) / (gamma + mu)
        beta <- (gamma + mu) * c(Z[-1], NA) /
          (S * c(Z[-1], NA))
      }
      list(S, I, beta)
    }
  )

  data
}

#' @rdname estimate-beta
#' @export
estimate_beta_si <- function(data, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "I0", "tgen")], envir = environment())

  # Preallocate memory for output
  data <- data[c("t", "Z", "B", "mu")]
  data[c("S", "I", "beta")] <- NA


  ## Estimate susceptible, infectious population sizes ...  -------------
  ## and transmission rate

  gamma <- 1 / tgen

  data[c("S", "I", "beta")] <- with(data,
    {
      S[1] <- S0
      I[1] <- I0
      for (i in 2:nrow(data)) {
        S[i] <- ((1 - 0.5 * mu[i-1] * 1) * S[i-1] + B[i] - Z[i]) /
          (1 + 0.5 * mu[i] * 1)
        I[i] <- ((1 - 0.5 * (gamma + mu[i-1]) * 1) * I[i-1] + Z[i]) /
          (1 + 0.5 * (gamma + mu[i]) * 1)
      }
      beta <- (Z + c(Z[-1], NA)) / (2 * S * I * 1)
      list(S, I, beta)
    }
  )

  structure(data, call = match.call(), arg_list = arg_list)
}

#' @rdname estimate-beta
#' @export
estimate_beta_sei <- function(data, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "E0", "I0", "tlat", "tinf")], envir = environment())

  # Preallocate memory for output
  data <- data[c("t", "Z", "B", "mu")]
  data[c("S", "E", "I", "beta")] <- NA


  ## Estimate susceptible, exposed, infectious population sizes ... ------
  ## and transmission rate

  sigma <- 1 / tlat
  gamma <- 1 / tinf

  data[c("S", "E", "I", "beta")] <- with(data,
    {
      S[1] <- S0
      E[1] <- E0
      I[1] <- I0
      for (i in 2:nrow(data)) {
        S[i] <- ((1-0.5*mu[i-1]*1)*S[i-1]+B[i]-Z[i]) /
          (1+0.5*mu[i]*1)
        E[i] <- ((1-0.5*(sigma+mu[i-1])*1)*E[i-1]+Z[i]) /
          (1+0.5*(sigma+mu[i])*1)
        I[i] <- ((1-0.5*(gamma+mu[i-1])*1)*I[i-1]+0.5*sigma*(E[i-1]+E[i])*1) /
          (1+0.5*(gamma+mu[i])*1)
      }
      beta <- (Z+c(Z[-1], NA)) / (2*S*I*1)
      list(S, E, I, beta)
    }
  )

  data
}
