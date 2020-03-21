#' Estimate
#' \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
#' using PTPI
#'
#' Using the method of peak-to-peak iteration (PTPI, see References),
#' `ptpi()` estimates the initial number of susceptibles
#' \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
#' from time series of incidence (*periodic*), births, and
#' natural mortality, observed at equally spaced time points
#' \ifelse{latex}{\out{$t_k = t_0 + k \Delta t$}}{\ifelse{html}{\out{<i>t<sub>k</sub></i> = <i>t</i><sub>0</sub>+<i>k&Delta;t</i>}}{t_k = t_0 + k*Dt}}
#' (for \ifelse{latex}{\out{$k = 0,\ldots,n$}}{\ifelse{html}{\out{<i>k</i> = 0,...,<i>n</i>}}{k = 0,...,n}}),
#' where
#' \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#' denotes the observation interval.
#'
#' @details
#' If `df$B` is undefined in the function call, then `df$B[i]` gets the
#' value `with(par_list, nu * hatN0 * 1)` for all `i`. If `df$mu` is
#' undefined in the function call, then `df$mu[i]` gets the value
#' `with(par_list, mu)` for all `i`.
#'
#' Missing values in `df[, c("Z", "B", "mu")]` are mostly not tolerated.
#' At the moment, `ptpi()` makes no effort to impute them, so imputation
#' must be done beforehand.
#'
#' @param df A data frame with numeric columns:
#'
#'   \describe{
#'     \item{`Z`}{Incidence. `Z[i]` is the number of infections between
#'       time
#'       \ifelse{latex}{\out{$t = t_{i-2}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;2</sub>}}{t = t_i-2}}
#'       and time
#'       \ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}.
#'       `Z` must be roughly periodic.
#'     }
#'     \item{`B`}{Births. `B[i]` is the number of births between time
#'       \ifelse{latex}{\out{$t = t_{i-2}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;2</sub>}}{t = t_i-2}}
#'       and time
#'       \ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}.
#'     }
#'     \item{`mu`}{Natural mortality rate. `mu[i]` is the rate at time
#'       \ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}
#'       expressed per unit
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
#'   }
#'
#'   `hatN0` and `nu` are optional if `B` is defined in `df`, and
#'   `mu` is optional if `mu` is defined in `df` (see Details).
#' @param a Integer scalar. Index of first peak in `df$Z`, possibly
#'   obtained via `get_peak_times()`.
#' @param b Integer scalar. Index of last peak in `df$Z` in phase
#'   with first peak, possibly obtained via `get_peak_times()`.
#' @param initial_S0_est Numeric scalar. An initial estimate of
#'   \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.
#' @param iter Integer scalar. The number of estimates of
#'   \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
#'   to generate before stopping.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{`S_mat`}{A numeric matrix with dimensions
#'     `c(nrow(df), iter+1)` containing the susceptible time series
#'     generated in each iteration.
#'   }
#'   \item{`S0`}{A numeric vector listing in order all `1+iter`
#'     estimates of
#'     \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.
#'     Equal to `S_mat[1, ]`.
#'   }
#'   \item{`S0_final`}{The final estimate of
#'     \ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.
#'     Equal to `S0[length(S0)]`.
#'   }
#'   \item{`SA`}{A numeric vector listing in order all `1+iter`
#'     estimates of
#'     \ifelse{latex}{\out{$S_\text{a} = S(t_\text{a})$}}{\ifelse{html}{\out{<i>S</i><sub>a</sub> = <i>S</i>(<i>t</i><sub>a</sub>)}}{S_a = S(t_a)}},
#'     where
#'     \ifelse{latex}{\out{$t_\text{a}$}}{\ifelse{html}{\out{<i>t</i><sub>a</sub>}}{t_a}}
#'     is the time point corresponding to row `a` in `df`.
#'     Equal to `S_mat[a, ]`.
#'   }
#'   \item{`SA_final`}{The final estimate of
#'     \ifelse{latex}{\out{$S_\text{a} = S(t_\text{a})$}}{\ifelse{html}{\out{<i>S</i><sub>a</sub> = <i>S</i>(<i>t</i><sub>a</sub>)}}{S_a = S(t_a)}}.
#'     Equal to `SA[length(SA)]`.
#'   }
#' }
#'
#' A list of the arguments of `ptpi()` is included as an attribute.
#'
#' @references
#' deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
#' time-varying transmission rates for infectious diseases.
#'
#' @md
#' @export
ptpi <- function(df = data.frame(), par_list = list(),
                 a, b,
                 initial_S0_est,
                 iter = 0L) {

## 1. Set-up -----------------------------------------------------------

# Assume constant vital rates if vital data were not supplied
if (is.null(df$B)) {
  df$B  <- with(par_list, nu * hatN0 * 1)
}
if (is.null(df$mu)) {
  df$mu <- with(par_list, mu)
}

# Preallocate memory for all susceptible time series,
# and initialize the first
S_mat <- matrix(NA, nrow = nrow(df), ncol = iter + 1)
S_mat[1, 1] <- initial_S0_est
S_mat[a, 1] <- initial_S0_est


## 2. Peak-to-peak iteration -------------------------------------------

for (j in seq_len(iter + 1)) {

  ## 2.(a) Update `SA` estimate

  # Reconstruct from index `a` to end
  for (i in (a+1):nrow(S_mat)) {
    S_mat[i, j] <- with(df[c("Z", "B", "mu")],
      {
        ((1 - 0.5 * mu[i-1] * 1) * S_mat[i-1,j] + B[i] - Z[i]) /
          (1 + 0.5 * mu[i] * 1)
      }
    )
  }
  if (j == iter + 1) {
    break
  }
  S_mat[a, j+1] <- S_mat[b, j]

  ## 2.(b) Update `S0` estimate

  # Reconstruct from index `a` to start (backwards in time)
  for (i in (a-1):1) {
    S_mat[i, j+1] <- with(df[c("Z", "B", "mu")],
      {
        ((1 + 0.5 * mu[i+1] * 1) * S_mat[i+1, j+1] - B[i+1] + Z[i+1]) /
          (1 - 0.5 * mu[i] * 1)
      }
    )
  }

}


out <- list(
  S_mat     = S_mat,
  S0       = S_mat[1, ],
  S0_final = S_mat[1, ncol(S_mat)],
  SA       = S_mat[a, ],
  SA_final = S_mat[a, ncol(S_mat)]
)
attr(out, "arg_list") <- as.list(environment())[names(formals(ptpi))]
out
}
