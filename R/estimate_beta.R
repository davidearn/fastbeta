#' \loadmathjax
#' Estimate time-varying transmission rates
#'
#' @description
#' Functions implementing the FC, S, SI, and SEI methods (see Algorithm)
#' for estimating time-varying transmission rates \mjseqn{\beta(t)} from
#' time series of incidence, births, and natural mortality with observations
#' at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}. The FC
#' and S methods are deprecated and outperformed by the more robust SI
#' and SEI methods. The SEI method accounts for a latent period and so
#' should be preferred over the SI method in practice. These 4 functions
#' should not be called directly, as they are all accessible through the
#' `method` argument of the more useful constructor function [fastbeta()],
#' which does checks on input and output and returns a fastbeta object
#' with associated methods.
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
#' ### SEI method
#' The susceptible, exposed, and infectious population sizes are estimated
#' recursively starting from the supplied initial values \mjseqn{S_0}
#' (`par_list$S0`), \mjseqn{E_0} (`par_list$E0`), and \mjseqn{I_0}
#' (`par_list$I0`):
#'
#' \mjsdeqn{\begin{align*} S_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} \mu_i \Delta t \big\rbrack S_i + B_{i+1} - Z_{i+1}}{1 + \frac{1}{2} \mu_{i+1} \Delta t}\,, \cr E_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\sigma + \mu_i) \Delta t \big\rbrack E_i + Z_i}{1 + \frac{1}{2} (\sigma + \mu_{i+1}) \Delta t}\,,\quad \sigma = 1 / t_\text{lat}\,, \cr I_{i+1} &= \frac{\big\lbrack 1 - \frac{1}{2} (\gamma + \mu_i) \Delta t \big\rbrack I_i + \frac{1}{2} \sigma (E_i + E_{i+1}) \Delta t}{1 + \frac{1}{2} (\gamma + \mu_{i+1}) \Delta t}\,,\quad \sigma = 1 / t_\text{lat}\,,\quad \gamma = 1 / t_\text{inf}\,. \end{align*}}
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
#' and SEI methods, strings of zeros in `df$Z` can lead to spurious
#' zeros and large numeric elements in column `beta` of the output.
#' Zeros in `df$Z` can be "imputed" (replaced with positive values
#' in a sensible way) to avoid these issues.
#'
#' Imputation is not carried out by `estimate_beta_fc()` and friends,
#' but *is* carried out by [fastbeta()].
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
#'     \item{`mu_i`}{\mjseqn{\lbrace\,\mu_i\,\rbrace}
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
#'       Number of exposed individuals at time \mjseqn{t = t_0}.
#'       SEI method only.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Number of infected (SI) or infectious (SEI) individuals
#'       at time \mjseqn{t = t_0}. SI and SEI methods only.
#'     }
#'   }
#'
#' @return
#' A data frame with numeric columns:
#'
#' \describe{
#'   \item{`t`}{\mjseqn{\lbrace\,t_i\,\rbrace} Time. Identical to `df$t`.}
#'   \item{`Z`}{\mjseqn{\lbrace\,Z_i\,\rbrace} Incidence. Identical to `df$Z`.}
#'   \item{`Z_agg`}{\mjseqn{\lbrace\,Z_i^\text{agg}\,\rbrace}
#'     Aggregated incidence. `Z_agg[i]` is the number of infections
#'     between times `t[i-round(tgen)+1]` and `t[i]`, for `i` in
#'     `seq(1 + round(tgen), nrow(df), by = round(tgen))`. `Z_agg[i]`
#'     is `NA` for all other `i`. FC method only.
#'   }
#'   \item{`B`}{\mjseqn{\lbrace\,B_i\,\rbrace} Births. Identical to `df$B`.}
#'   \item{`B_agg`}{\mjseqn{\lbrace\,B_i^\text{agg}\,\rbrace}
#'     Aggregated births. `B_agg[i]` is the number of births
#'     between times `t[i-round(tgen)+1]` and `t[i]`, for `i` in
#'     `seq(1 + round(tgen), nrow(df), by = round(tgen))`. `B_agg[i]`
#'     is `NA` for all other `i`. FC method only.
#'   }
#'   \item{`mu`}{\mjseqn{\lbrace\,\mu_i\,\rbrace}
#'     Per capita natural mortality rate. Identical to `df$mu`.
#'     S, SI, and SEI methods only.
#'   }
#'   \item{`S`}{\mjseqn{\lbrace\,S_i\,\rbrace}
#'     Estimated number of susceptible individuals.
#'   }
#'   \item{`E`}{\mjseqn{\lbrace\,E_i\,\rbrace}
#'     Estimated number of exposed individuals. SEI method only.
#'   }
#'   \item{`I`}{\mjseqn{\lbrace\,I_i\,\rbrace}
#'     Estimated number of infected (FC, S, SI) or infectious (SEI)
#'     individuals.
#'   }
#'   \item{`beta`}{\mjseqn{\lbrace\,\beta_i\,\rbrace}
#'     Estimated transmission rate, expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infected (FC, S, SI)
#'     or infectious (SEI) individual.
#'   }
#' }
#'
#' The data frame has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(estimate_beta_x, arg_list)`.
#'
#' @examples
#' # Simulate time series data using an SIR model
#' # to test the FC, S, and SI methods
#' par_list_sir <- make_par_list(model = "sir")
#' df_sir <- make_data(par_list_sir, with_ds = TRUE, model = "sir")
#' head(df_sir)
#'
#' # Simulate time series data using an SEIR model
#' # to test the SEI method
#' par_list_seir <- make_par_list(model = "seir")
#' df_seir <- make_data(par_list_seir, with_ds = TRUE, model = "seir")
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
#' eb_out <- mapply(
#'   function(f, df, par_list) f(df, par_list),
#'   f = f_list,
#'   df = list(df_sir, df_sir, df_sir, df_seir),
#'   par_list = list(
#'     within(par_list_sir, tgen <- tlat + tinf), #' fc
#'     within(par_list_sir, tgen <- tlat + tinf), #' s
#'     within(par_list_sir, tgen <- tlat + tinf), #' si
#'     par_list_seir                              #' sei
#'   ),
#'   SIMPLIFY = FALSE
#' )
#' lapply(eb_out, head)
#'
#' op <- par(mfrow = c(3, 1), mar = c(2, 5, 1, 1), oma = c(3, 0, 2, 0))
#'
#' # FC method fails rapidly as a result of ignoring
#' # susceptible mortality
#' plot(S ~ t, df_sir, type = "l", lwd = 4, ylim = c(43, 83) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, na.omit(eb_out$fc), lwd = 2, col = "seagreen")
#' title(main = "FC method", line = 1, xpd = NA)
#' plot(I ~ t, df_sir, type = "l", lwd = 4, ylim = c(0, 3.5) * 1e03,
#'      xlab = "", ylab = "Infected")
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
#'      xlab = "", ylab = "Infected")
#' lines(I ~ t, eb_out$s, lwd = 2, col = "mediumvioletred")
#' plot(beta ~ t, df_sir, type = "l", lwd = 4, ylim = c(1.0, 1.4) * 1e-05,
#'      xlab = "", ylab = "Transmission rate")
#' lines(beta ~ t, eb_out$s, lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' # SI method does well
#' plot(S ~ t, df_sir, type = "l", lwd = 4, ylim = c(43, 58) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, eb_out$si, lwd = 2, col = "seagreen")
#' title(main = "SI method", line = 1, xpd = NA)
#' plot(I ~ t, df_sir, type = "l", lwd = 4, ylim = c(0, 3.5) * 1e03,
#'      xlab = "", ylab = "Infected")
#' lines(I ~ t, eb_out$si, lwd = 2, col = "mediumvioletred")
#' plot(beta ~ t, df_sir, type = "l", lwd = 4, ylim = c(1.0, 1.4) * 1e-05,
#'      xlab = "", ylab = "Transmission rate")
#' lines(beta ~ t, eb_out$si, lwd = 2, col = "slateblue")
#' title(xlab = "Time (units dt)", line = 3, xpd = NA)
#'
#' # SEI method does well, but see note below
#' plot(S ~ t, df_seir, type = "l", lwd = 4, ylim = c(43, 58) * 1e03,
#'      xlab = "", ylab = "Susceptible")
#' lines(S ~ t, eb_out$sei, lwd = 2, col = "seagreen")
#' title(main = "SEI method", line = 1, xpd = NA)
#' plot(I ~ t, df_seir, type = "l", lwd = 4, ylim = c(0, 2) * 1e03,
#'      xlab = "", ylab = "Infectious")
#' lines(I ~ t, eb_out$sei, lwd = 2, col = "mediumvioletred")
#' plot(beta ~ t, df_seir, type = "l", lwd = 4, ylim = c(1.7, 2.3) * 1e-05,
#'      xlab = "", ylab = "Transmission rate")
#' lines(beta ~ t, eb_out$sei, lwd = 2, col = "slateblue")
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
estimate_beta_fc <- function(df, par_list) {
  ## Set-up --------------------------------------------------------------

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


  ## Aggregate incidence, births -----------------------------------------

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


  ## Estimate susceptibles, infecteds, transmission rate -----------------

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

#' @rdname estimate-beta
#' @export
estimate_beta_s <- function(df, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "tgen")], envir = environment())

  # Preallocate memory for output
  df <- df[c("t", "Z", "B", "mu")]
  df[c("S", "I", "beta")] <- NA


  ## Estimate susceptibles, infecteds, transmission rate -----------------

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

#' @rdname estimate-beta
#' @export
estimate_beta_si <- function(df, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "I0", "tgen")], envir = environment())

  # Preallocate memory for output
  df <- df[c("t", "Z", "B", "mu")]
  df[c("S", "I", "beta")] <- NA


  ## Estimate susceptibles, infecteds, transmission rate -----------------

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

#' @rdname estimate-beta
#' @export
estimate_beta_sei <- function(df, par_list) {
  ## Set-up --------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(par_list[c("S0", "E0", "I0", "tlat", "tinf")], envir = environment())

  # Preallocate memory for output
  df <- df[c("t", "Z", "B", "mu")]
  df[c("S", "E", "I", "beta")] <- NA


  ## Estimate susceptible, exposed, infectious population sizes ... ------
  ## and transmission rate

  sigma <- 1 / tlat
  gamma <- 1 / tinf

  df[c("S", "E", "I", "beta")] <- with(df,
    {
      S[1] <- S0
      E[1] <- E0
      I[1] <- I0
      for (i in 2:nrow(df)) {
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

  structure(df, call = match.call(), arg_list = arg_list)
}
