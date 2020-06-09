#' Estimate time-varying transmission rates (S method)
#'
#' `estimate_beta_S()` applies the S method (see References)
#' to estimate the time-varying transmission rate
#' \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#' from time series of reported incidence, births, and
#' natural mortality, observed at equally spaced time points
#' \ifelse{latex}{\out{$t_k = t_0 + k \Delta t$}}{\ifelse{html}{\out{<i>t<sub>k</sub></i> = <i>t</i><sub>0</sub>+<i>k&Delta;t</i>}}{t_k = t_0 + k*Dt}}
#' (for \ifelse{latex}{\out{$k = 0,\ldots,n$}}{\ifelse{html}{\out{<i>k</i> = 0,...,<i>n</i>}}{k = 0,...,n}}),
#' where
#' \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#' denotes the observation interval.
#'
#' @section Mock vital data:
#' If `df$B` is undefined in the function call, then `df$B[i]`
#' gets the value `with(par_list, nu * hatN0 * 1)` for all `i`.
#' If `df$mu` is undefined the function call, then `df$mu[i]`
#' gets the value `with(par_list, mu)` for all `i`.
#'
#' @section Missing data:
#' Missing values in `df[, c("C", "B", "mu")]` are not tolerated
#' by the S method. They are imputed via linear interpolation
#' between observed values. If there are no observations before
#' the first missing value, then complete imputation is impossible.
#' In this case, the S method may fail: columns `S` and `beta` in
#' the output may be filled with `NA`.
#'
#' Zeros in `df$C` cause divide-by-zero errors. To prevent these
#' errors, zeros are imputed like missing values. If there are
#' no nonzero observations before the first zero, then complete
#' imputation is impossible. In this case, the S method may fail,
#' but only locally: column `beta` in the output may contain some
#' `NaN` and `Inf`, but give estimates everywhere else.
#'
#' @param df A data frame with numeric columns:
#'
#'   \describe{
#'     \item{`t`}{Time. `t[i]` is equal to
#'       \ifelse{latex}{\out{$t_{i-1} = t_0 + (i-1) \Delta t$}}{\ifelse{html}{\out{<i>t</i><sub><i>i</i>&minus;1</sub> = <i>t</i><sub>0</sub>+(<i>i</i>&minus;1)<i>&Delta;t</i>}}{t_i-1 = t_0 + (i-1)*Dt}}
#'       in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}},
#'       so that `t[i] - t[i-1]` is equal to 1.
#'     }
#'     \item{`C`}{Reported incidence. `C[i]` is the number of cases
#'       reported between times `t[i-1]` and `t[i]`.
#'     }
#'     \item{`B`}{Births. `B[i]` is the number of births between times
#'       `t[i-1]` and `t[i]`.
#'     }
#'     \item{`mu`}{Natural mortality rate. `mu[i]` is the rate at time
#'       `t[i]` expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       and per capita.
#'     }
#'   }
#'
#'   `B` is optional if `hatN0` and `nu` are defined in `par_list`, and
#'   `mu` is optional if `mu` is defined in `par_list` (see Details).
#' @param par_list A list containing:
#'
#'   \describe{
#'     \item{`prep`}{\[ \ifelse{latex}{\out{$p_\text{rep}$}}{\ifelse{html}{\out{<i>p</i><sub>rep</sub>}}{p_rep}} \]
#'       Case reporting probability.
#'     }
#'     \item{`trep`}{\[ \ifelse{latex}{\out{$t_\text{rep}$}}{\ifelse{html}{\out{<i>t</i><sub>rep</sub>}}{t_rep}} \]
#'       Case reporting delay in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     }
#'     \item{`S0`}{\[ \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}} \]
#'       Number of susceptibles at time
#'       \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'     }
#'     \item{`hatN0`}{\[ \ifelse{latex}{\out{$\widehat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}} \]
#'       Population size at time
#'       \ifelse{latex}{\out{$t = 0$}}{\ifelse{html}{\out{<i>t</i> = 0}}{t = 0}}
#'       years.
#'     }
#'     \item{`nu`}{\[ \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}} \]
#'       Birth rate expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       and relative to
#'       \ifelse{latex}{\out{$\hat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}}
#'       (if modeled as constant).
#'     }
#'     \item{`mu`}{\[ \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}} \]
#'       Natural mortality rate expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       and per capita (if modeled as constant).
#'     }
#'     \item{`tgen`}{\[ \ifelse{latex}{\out{$t_\text{gen}$}}{\ifelse{html}{\out{<i>t</i><sub>gen</sub>}}{t_gen}} \]
#'       Mean generation interval of the disease of interest in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     }
#'   }
#'
#'   `hatN0` and `nu` are optional if `B` is defined in `df`, and
#'   `mu` is optional if `mu` is defined in `df` (see Details).
#'
#' @return
#' A data frame with numeric columns:
#'
#' \describe{
#'   \item{`t`}{Time. Identical to `df$t`.}
#'   \item{`C`}{Reported incidence, imputed. Identical to `df$C`,
#'     except with missing values and zeros imputed (see Details).
#'   }
#'   \item{`Z`}{Incidence. `Z[i]` is the estimated
#'     number of infections between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`B`}{Births, imputed. Identical to `df$B` (if supplied),
#'     except with missing values imputed (see Details).
#'   }
#'   \item{`mu`}{Natural mortality rate, imputed. Identical to `df$mu`
#'     (if supplied), except with missing values imputed (see Details).
#'   }
#'   \item{`S`}{Number of susceptibles. `S[i]` is the estimated
#'     number of susceptibles at time `t[i]`.
#'   }
#'   \item{`I`}{Number of infecteds. `I[i]` is the estimated
#'     number of infecteds at time `t[i]`.
#'   }
#'   \item{`beta`}{Transmission rate. `beta[i]` is the estimated
#'     transmission rate at time `t[i]` expressed per unit
#'     \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'     per susceptible per infected.
#'   }
#' }
#'
#' It possesses `par_list` as an attribute.
#'
#' @examples
#' # Simulate a reported incidence time series using
#' # a seasonally forced transmission rate
#' par_list <- make_par_list(dt_weeks = 1, epsilon = 0.5, prep = 0.5)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7, # 20 years is ~1042 weeks
#'   with_dem_stoch = TRUE,
#'   seed = 5
#' )
#' head(df)
#'
#' # Estimate incidence, susceptibles, infecteds,
#' # and the seasonally forced transmission rate
#' df_S <- estimate_beta_S(df, par_list)
#' head(df_S)
#'
#' # Fit a smooth loess curve to the transmission rate
#' # time series
#' loess_fit <- loess(
#'   formula   = beta ~ t,
#'   data      = df_S,
#'   span      = 65 / nrow(df_S),
#'   degree    = 2,
#'   na.action = "na.exclude"
#' )
#' df_S$beta_loess <- predict(loess_fit)
#'
#' # Inspect
#' df_S$t_years <- df$t_years
#' plot(S ~ t_years, df, type = "l", ylim = c(43, 58) * 1e03)
#' lines(S ~ t_years, df_S, col = "red")
#' plot(beta ~ t_years, df, type = "l", ylim = c(0.95, 1.25) * 1e-05)
#' lines(beta_loess ~ t_years, df_S, col = "red")
#' 
#' @references
#' deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
#' time-varying transmission rates for infectious diseases.
#'
#' @md
#' @export
estimate_beta_S <- function(df       = data.frame(),
                            par_list = list()) {

## 1. Set-up -----------------------------------------------------------

# Load necessary elements of `par_list` into the execution environment
list2env(
  par_list[c("prep", "trep", "S0", "hatN0", "nu", "mu", "tgen")],
  envir = environment()
)

# Logical: Is `df` missing a column `B`? What about `mu`?
B_was_not_def <- is.null(df$B)
mu_was_not_def <- is.null(df$mu)

# Preallocate memory for output. Assume constant vital rates,
# if necessary.
df <- data.frame(
  t     = df$t,
  C     = df$C,
  Z     = NA,
  B     = if (B_was_not_def) nu * hatN0 * 1 else df$B,
  mu    = if (mu_was_not_def) mu else df$mu,
  S     = NA,
  I     = NA,
  beta  = NA
)


## 2. Impute missing values in reported incidence, ... -----------------
##    births, and natural mortality rate

df[c("C", "B", "mu")] <- lapply(df[c("C", "B", "mu")],
  function(x) {
    if (!any(is.na(x))) {
      return(x)
    }

    # Indices of missing values
    ind_na <- which(is.na(x))

    # Function that performs linear interpolation
    # between time points where `x` is observed
    impute_na <- stats::approxfun(
      x      = df$t[-ind_na],
      y      = x[-ind_na],
      method = "linear",
      rule   = 1 # return `NA` outside range of (argument) `x`
    )

    # Impute missing values. Missing values outside of
    # the range of observed data are retained as `NA`.
    replace(x, ind_na, impute_na(df$t[ind_na]))
  }
)


## 3. Impute zeros in reported incidence -------------------------------

df$C <- local(
  {
    if (!any(df$C == 0, na.rm = TRUE)) {
      return(df$C)
    }

    # Indices of positive values
    ind_pos <- which(df$C > 0)

    # Indices of zeros between two positive values
    ind_zero <- which(df$C == 0)
    ind_zero <- ind_zero[
      ind_zero > min(ind_pos) &
      ind_zero < max(ind_pos)
    ]

    # Function that performs linear interpolation
    # between time points where `df$C` is positive
    impute_zero <- stats::approxfun(
      x      = df$t[ind_pos],
      y      = df$C[ind_pos],
      method = "linear",
      rule   = 1 # return `NA` outside range of `x`
    )

    # Impute zeros. Zeros outside of the range of
    # positive data are retained as zeros due to
    # the definition of `ind_zero`.
    replace(df$C, ind_zero, impute_zero(df$t[ind_zero]))
  }
)


## NOTE: In the best case, there are no longer missing values or zeros
##       in `df$C`. In the worst case, `df$C` now looks like
##       `c(NA,...,NA,0,...,0,+,...,+,0,...,0,NA,...,NA)`,
##       where `+` denotes a positive number.


## 4. Estimate incidence -----------------------------------------------

trepr <- round(trep)
df$Z <- c(
  (1 / prep) * df$C[(trepr+1):nrow(df)],
  rep(NA, trepr)
)


## 5. Estimate susceptibles, infecteds, transmission rate --------------

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


## NOTE: If too many `NA` were retained at the start of `df$C`,
##       `df$B`, or `df$mu`, then `df$S` and `df$beta` will be
##       filled with `NA`. If zeros were retained in `df$C`,
##       then `NaN` or `Inf` will appear in `df$beta` whenever
##       a divide-by-zero error occurs.


## 6. Warn if `S` is ever negative -------------------------------------

if (any(df$S < 0, na.rm = TRUE)) {
  # May have underestimated `prep` or `nu`
  new_par_vals <- paste0(
    "\n* `prep` > ", sprintf("%.3f", prep),
    if (B_was_not_def) paste0("\n* `nu` > ", sprintf("%.3e", nu))
  )
  warning(
    "S method: `S[i]` < 0 for some `i`. Retry with:",
    new_par_vals,
    call. = FALSE
  )
}


## 7. Warn if `beta` is ever `NaN` or `Inf` ----------------------------

if (any(is.nan(df$beta) | is.infinite(df$beta))) {
  warning(
    "S method: `beta[i]` is `NaN` or `Inf` for some `i`. \n",
    "Is the first or last observation in `df$C` a zero?",
    call. = FALSE
  )
}


attr(df, "par_list") <- par_list
df
}
