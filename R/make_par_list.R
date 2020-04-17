#' Create a list of parameter values
#'
#' @description
#' `make_par_list()` creates a list of parameter values that may be
#' passed as an argument of:
#'
#' * [make_data()] to simulate time series data.
#' * [estimate_beta_FC()], [estimate_beta_S()], and [estimate_beta_SI()]
#'   to estimate time-varying transmission rates
#'   \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}.
#' * [ptpi()] to estimate the initial number of susceptibles
#'   \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}}.
#'
#' As these methods deal with equally spaced time series data
#' with fixed observation interval
#' \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}},
#' `make_par_list()` defines time (or rate) parameters in units
#' (or per unit)
#' \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'
#' @details
#' `make_par_list()` enforces the identity
#'
#' \ifelse{latex}{\out{$\mathcal{R}_0 = \frac{\nu_\text{c} \widehat{N}_0}{\mu_\text{c}} \cdot \frac{\langle\beta\rangle}{\gamma + \mu}$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub> = (<i>&nu;</i><sub>c</sub> <i>&Ntilde;</i> / <i>&mu;</i><sub>c</sub>)(&langle;<i>&beta;</i>&rangle; / (<i>&gamma;</i> + <i>&mu;</i><sub>c</sub>))}}{calR_0 = ((nu_c*hatN0)/mu_c)*(<beta>/(gamma + mu))}}
#'
#' as follows. If exactly one of `Rnaught` and `beta_mean` is defined
#' in the function call, then the parameter not defined is internally
#' assigned the value satisfying the identity. If both are defined in
#' the function call, then the value of `beta_mean` is replaced with the
#' value satisfying the identity.
#'
#' If `N0`, `S0`, or `I0` are not defined in the function call, then
#' `make_par_list()` numerically integrates the system of SIR equations
#'
#' \ifelse{latex}{
#'   \out{
#'     \begin{array}[rlc]
#'       S' & = & \nu_\text{c} \widehat{N}_0 - \beta(t) S I - \mu_\text{c} S \\
#'       I' & = & \beta(t) S I - \gamma I - \mu_\text{c} I \\
#'       R' & = & \gamma I - \mu_\text{c} R
#'     \end{array}
#'   }
#' }{
#'   \ifelse{html}{
#'     \out{
#'       <i>S</i>&prime; = <i>&nu;</i><sub>c</sub><i>&Ntilde;</i><sub>0</sub> &minus; <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&mu;</i><sub>c</sub><i>S</i><br>
#'       <i>I</i>&prime; = <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&gamma;I</i> &minus; <i>&mu;</i><sub>c</sub><i>I</i><br>
#'       <i>R</i>&prime; = <i>&gamma;I</i> &minus; <i>&mu;</i><sub>c</sub><i>R</i>
#'     }
#'   }{
#'     S' = nu_c*hatN_0 - beta(t)*S*I - mu_c*S \cr
#'     I' = beta(t)*S*I - gamma*I - mu_c*I \cr
#'     R' = gamma*I - mu_c*R
#'   }
#' }
#'
#' with
#' \ifelse{latex}{\out{$\gamma = 1 / t_\text{gen}$}}{\ifelse{html}{\out{<i>&gamma;</i> = 1 / <i>t</i><sub>gen</sub>}}{gamma = 1/t_gen}}
#' and
#' \ifelse{latex}{\out{$\beta(t) = \langle\beta\rangle \left(1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}}\right)\right)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>) = &langle;<i>&beta;</i>&rangle; (1 + <i>&alpha;</i> cos(2<i>&pi;t</i> / (1 year)))}}{beta(t) = <beta>*(1 + alpha*cos(2*pi*t/(1 year)))}}
#' between times
#' \ifelse{latex}{\out{$t = 0$}}{\ifelse{html}{\out{<i>t</i> = 0}}{t = 0}}
#' years and
#' \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}},
#' taking for the initial state
#' \ifelse{latex}{\out{$(S(0),I(0),R(0))$}}{\ifelse{html}{\out{(<i>S</i>(0),<i>I</i>(0),<i>R</i>(0))}}{(S(0),I(0),R(0))}}
#' the endemic equilibrium of the system with
#' \ifelse{latex}{\out{$\beta(t) \equiv \langle\beta\rangle$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>) &equiv; &langle;<i>&beta;</i>&rangle;}}{beta(t) = <beta>}}
#' and
#' \ifelse{latex}{\out{$\nu_\text{c} = \mu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;</i><sub>c</sub> = <i>&mu;</i><sub>c</sub>}}{nu_c = mu_c}}.
#' Then `make_par_list()` defines `N0`, `S0`, and `I0` as follows:
#'
#' \describe{
#'   \item{`N0`}{The value of
#'     \ifelse{latex}{\out{$S(t_0)+I(t_0)+R(t_0)$}}{\ifelse{html}{\out{<i>S</i>(<i>t</i><sub>0</sub>)+<i>I</i>(<i>t</i><sub>0</sub>)+<i>R</i>(<i>t</i><sub>0</sub>)}}{S(t_0)+I(t_0)+R(t_0)}}.
#'   }
#'   \item{`S0`}{The value of
#'     \ifelse{latex}{\out{$S(t_0)$}}{\ifelse{html}{\out{<i>S</i>(<i>t</i><sub>0</sub>)}}{S(t_0)}}.
#'   }
#'   \item{`I0`}{The value of
#'     \ifelse{latex}{\out{$I(t_0)$}}{\ifelse{html}{\out{<i>I</i>(<i>t</i><sub>0</sub>)}}{I(t_0)}}.
#'   }
#' }
#'
#' @param dt_weeks \[ \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}} \]
#'   Observation interval in weeks.
#' @param t0 \[ \ifelse{latex}{\out{$t_0$}}{\ifelse{html}{\out{<i>t</i><sub>0</sub>}}{t_0}} \]
#'   Time of the first observation in units
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#' @param prep \[ \ifelse{latex}{\out{$p_\text{rep}$}}{\ifelse{html}{\out{<i>p</i><sub>rep</sub>}}{p_rep}} \]
#'   Case reporting probability.
#' @param trep \[ \ifelse{latex}{\out{$t_\text{rep}$}}{\ifelse{html}{\out{<i>t</i><sub>rep</sub>}}{t_rep}} \]
#'   Case reporting delay in units
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#' @param hatN0 \[ \ifelse{latex}{\out{$\widehat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}} \]
#'   Population size at time
#'   \ifelse{latex}{\out{$t = 0$}}{\ifelse{html}{\out{<i>t</i> = 0}}{t = 0}}
#'   years (see Details).
#' @param N0 \[ \ifelse{latex}{\out{$N_0$}}{\ifelse{html}{\out{<i>N</i><sub>0</sub>}}{N_0}} \]
#'   Population size at time
#'   \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'   Can be set to `NULL` (see Details).
#' @param S0 \[ \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}} \]
#'   Number of susceptibles at time
#'   \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'   Can be set to `NULL` (see Details).
#' @param I0 \[ \ifelse{latex}{\out{$I_0$}}{\ifelse{html}{\out{<i>I</i><sub>0</sub>}}{I_0}} \]
#'   Number of infecteds at time
#'   \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'   Can be set to `NULL` (see Details).
#' @param nu \[ \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}} \]
#'   Birth rate expressed per unit
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'   and relative to
#'   \ifelse{latex}{\out{$\hat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}}
#'   (if modeled as constant).
#' @param mu \[ \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}} \]
#'   Natural mortality rate expressed per unit
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'   and per capita (if modeled as constant).
#' @param tgen \[ \ifelse{latex}{\out{$t_\text{gen}$}}{\ifelse{html}{\out{<i>t</i><sub>gen</sub>}}{t_gen}} \]
#'   Mean generation interval of the disease of interest in units
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#' @param Rnaught \[ \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}} \]
#'   Basic reproduction number of the disease of interest. Should be set
#'   to `NULL` when specifying `beta_mean` (see Details).
#' @param beta_mean \[ \ifelse{latex}{\out{$\langle\beta\rangle$}}{\ifelse{html}{\out{&langle;<i>&beta;</i>&rangle;}}{<beta>}} \]
#'   Mean of the seasonally forced transmission rate
#'   \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#'   expressed per unit
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'   per susceptible per infected. Should be set to `NULL` when
#'   specifying `Rnaught` (see Details).
#' @param alpha \[ \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}} \]
#'   Amplitude of the seasonally forced transmission rate
#'   \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#'   relative to the mean.
#' @param epsilon \[ \ifelse{latex}{\out{$\epsilon$}}{\ifelse{html}{\out{<i>&straightepsilon;</i>}}{epsilon}} \]
#'   Standard deviation of the random phase shift in the seasonally
#'   forced transmission rate
#'   \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}.
#'
#' @return
#' A list of the arguments of `make_par_list()`, including values for
#' `Rnaught`, `beta_mean`, `N0`, `S0`, and `I0` if not defined in the
#' function call (see Details).
#'
#' @examples
#' # All arguments have default values
#' par_list <- make_par_list()
#' unlist(par_list)
#'
#' # Time and rate parameters must be specified
#' # in terms of the observation interval
#' dt_weeks <- 1
#' nu_peryear <- 0.04
#' tgen_days <- 13
#' par_list <- make_par_list(
#'   dt_weeks = dt_weeks,
#'   nu       = nu_peryear * (7 / 365) * dt_weeks,
#'   tgen     = tgen_days * (1 / 7) / dt_weeks
#' )
#' unlist(par_list)
#' 
#' @md
#' @export
make_par_list <- function(dt_weeks  = 1,
                          t0        = 2000 * (365 / 7) / dt_weeks,
                          prep      = 1,
                          trep      = 0,
                          hatN0     = 1e06,
                          N0        = NULL,
                          S0        = NULL,
                          I0        = NULL,
                          nu        = 0.04 * (7 / 365) * dt_weeks,
                          mu        = 0.04 * (7 / 365) * dt_weeks,
                          tgen      = 13 * (1 / 7) / dt_weeks,
                          Rnaught   = 20,
                          beta_mean = NULL,
                          alpha     = 0.08,
                          epsilon   = 0.5) {

# Derived quantities
gamma <- 1 / tgen
one_year <- (365 / 7) / dt_weeks

# If `Rnaught` is defined but not `beta_mean`
if (!is.null(Rnaught) && is.null(beta_mean)) {
  beta_mean <- (mu / (nu * hatN0)) * Rnaught * (gamma + mu)
# If `beta_mean` is defined but not Rnaught`
} else if (!is.null(beta_mean) && is.null(Rnaught)) {
  Rnaught <- ((nu * hatN0) / mu) * beta_mean / (gamma + mu)
# If `Rnaught` and `beta_mean` are both defined
} else if (!is.null(Rnaught) && !is.null(beta_mean)) {
  beta_mean <- (mu / (nu * hatN0)) * Rnaught * (gamma + mu)
# Otherwise
} else {
  stop(
    "At least one of `Rnaught` and `beta_mean` must be specified.",
    call. = FALSE
  )
}


# If `N0`, `S0`, or `I0` is not defined, then
# obtain a value from the state of a system of
# SIR equations after a transient of length `t0`
if (is.null(N0) || is.null(S0) || is.null(I0)) {

  # Initial state
  x_init <- c(
    S = hatN0 * (1 / Rnaught),
    I = hatN0 * (1 - 1 / Rnaught) * (mu / (gamma + mu)),
    R = hatN0 * (1 - 1 / Rnaught) * (gamma / (gamma + mu))
  )

  # Time points
  t_out <- 0:floor(t0)

  # Seasonally forced transmission rate
  beta <- function(t) {
    beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
  }

  # System of SIR equations
  compute_sir_rates <- function(t, y, parms) {
    with(as.list(c(y, parms)),
      {
        dS <- nu * hatN0 - beta(t) * S * I - mu * S
        dI <- beta(t) * S * I - gamma * I - mu * I
        dR <- gamma * I - mu * R
        list(c(dS, dI, dR))
      }
    )
  }

  # Numerical solution of system
  tsir_mat <- deSolve::lsoda(
    x_init, t_out, compute_sir_rates,
    parms    = NULL,
    rtol     = 1e-16,
    maxsteps = 1e05
  )

  # Assign final values of `S+I+R`, `S`, and `I`, if necessary
  if (is.null(N0)) {
    N0 <- sum(tsir_mat[nrow(tsir_mat), c("S", "I", "R")])
  }
  if (is.null(S0)) {
    S0 <- as.numeric(tsir_mat[nrow(tsir_mat), "S"])
  }
  if (is.null(I0)) {
    I0 <- as.numeric(tsir_mat[nrow(tsir_mat), "I"])
  }

}


as.list(environment())[names(formals(make_par_list))]
}
