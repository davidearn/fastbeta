#' \loadmathjax
#' Create a list of parameter values for simulations
#'
#' @description
#' `make_par_list()` creates a list of parameter values to be passed as
#' an argument of [make_data()], which simulates epidemic time series data
#' using an SIR model.
#'
#' @param dt_weeks \mjseqn{\lbrack \Delta t \rbrack}
#'   Numeric scalar. Observation interval in weeks.
#' @param t0 \mjseqn{\lbrack t_0 \rbrack}
#'   Numeric scalar. Time of the first observation
#'   in units \mjseqn{\Delta t}.
#' @param prep \mjseqn{\lbrack p_\text{rep} \rbrack}
#'   Numeric scalar. Case reporting probability.
#' @param trep \mjseqn{\lbrack m \rbrack}
#'   Numeric scalar. Case reporting delay in units \mjseqn{\Delta t}.
#' @param hatN0 \mjseqn{\lbrack \widehat{N}_0 \rbrack}
#'   Numeric scalar. Population size at time \mjseqn{t = 0} years.
#' @param N0 \mjseqn{\lbrack N_0 \rbrack}
#'   Numeric scalar. Population size at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details 2).
#' @param S0 \mjseqn{\lbrack S_0 \rbrack}
#'   Numeric scalar. Number of susceptibles at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details 2).
#' @param I0 \mjseqn{\lbrack I_0 \rbrack}
#'   Numeric scalar. Number of infecteds at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details 2).
#' @param nu \mjseqn{\lbrack \nu \rbrack}
#'   Numeric scalar. Birth rate expressed per unit \mjseqn{\Delta t}
#'   and relative to \mjseqn{\widehat{N}_0}.
#' @param mu \mjseqn{\lbrack \mu \rbrack}
#'   Numeric scalar. Natural mortality rate expressed per unit
#'   \mjseqn{\Delta t} and per capita.
#' @param tgen \mjseqn{\lbrack t_\text{gen} \rbrack}
#'   Numeric scalar. Mean generation interval of the disease
#'   of interest in units \mjseqn{\Delta t}.
#' @param Rnaught \mjseqn{\lbrack \mathcal{R}_0 \rbrack}
#'   Numeric scalar. Basic reproduction number of the disease
#'   of interest. If specified (not `NULL`), then `beta_mean` is
#'   calculated internally as a function of `Rnaught` (see Details 1).
#' @param beta_mean \mjseqn{\lbrack \langle\beta\rangle \rbrack}
#'   Numeric scalar. Mean (long-term average) of the seasonally
#'   forced transmission rate \mjseqn{\beta(t)} expressed per unit
#'   \mjseqn{\Delta t} per susceptible per infected.
#'   Ignored if `Rnaught` is specified (not `NULL`) (see Details 1).
#' @param alpha \mjseqn{\lbrack \alpha \rbrack}
#'   Numeric scalar. Amplitude of the seasonally forced transmission
#'   rate \mjseqn{\beta(t)} relative to the mean.
#' @param epsilon \mjseqn{\lbrack \epsilon \rbrack}
#'   Numeric scalar. Standard deviation of the random phase shift
#'   introduced to the seasonally forced transmission rate
#'   \mjseqn{\beta(t)}.
#' @param ode_control A list of optional arguments of [deSolve::ode()],
#'   specifying options for numerical integration, such as `method`,
#'   `rtol`, and `atol` (see Details 2).
#'
#' @return
#' A list of the arguments of `make_par_list()` (excluding
#' `ode_control`) with values for `Rnaught`, `beta_mean`,
#' `N0`, `S0`, and `I0` if not supplied in the function call
#' (see Details 1 and 2).
#'
#' @details
#' # Details
#'
#' ## 1. Concordance of \mjseqn{\mathcal{R}_0} and \mjseqn{\langle\beta\rangle}
#' 
#' `make_par_list()` enforces the identity
#'
#' \mjsdeqn{\mathcal{R}_0 = \frac{\nu \widehat{N}_0}{\mu} \cdot \frac{\langle\beta\rangle}{\gamma + \mu}\,,}
#'
#' where \mjseqn{\gamma = 1 / t_\text{gen}}. This means that it calculates
#' \mjseqn{\langle\beta\rangle} as a function of \mjseqn{\mathcal{R}_0}, unless
#' the user explicitly sets `Rnaught = NULL`, in which case it calculates
#' \mjseqn{\mathcal{R}_0} as a function of \mjseqn{\langle\beta\rangle}. Hence
#' argument `beta_mean` is optional if `Rnaught` is specified, but mandatory
#' otherwise.
#' 
#' ## 2. Missing \mjseqn{N_0}, \mjseqn{S_0}, or \mjseqn{I_0}
#'
#' If \mjseqn{N_0}, \mjseqn{S_0}, or \mjseqn{I_0} (`N0`, `S0`, or `I0`)
#' is not specified in the function call, then, with a call to
#' [deSolve::ode()], `make_par_list()` numerically integrates the
#' system of SIR equations
#'
#' \mjsdeqn{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu \widehat{N}_0 - \beta(t) S I - \mu S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu I \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu R \end{align}}
#' 
#' with \mjseqn{\gamma = 1 / t_\text{gen}} and
#' 
#' \mjsdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}}\right) \right\rbrack}
#' 
#' between times \mjseqn{t = 0} years and \mjseqn{t = t_0}, taking
#'
#' \mjsdeqn{\begin{bmatrix} S(0) \cr I(0) \cr R(0) \end{bmatrix} = \widehat{N}_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu}{\gamma + \mu} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\gamma}{\gamma + \mu} \end{bmatrix}}
#'
#' for the initial state. This is the endemic equilibrium of the above system
#' with constant transmission rate \mjseqn{\beta(t) \equiv \langle\beta\rangle},
#' scaled by \mjseqn{\mu/\nu} to enforce
#' \mjseqn{S(0) + I(0) + R(0) = \widehat{N}_0}
#'
#' as the initial population size. Then `make_par_list()` assigns `N0`,
#' `S0`, and `I0` (only those not specified in the function call) the
#' values
#'
#' \describe{
#'   \item{`N0`}{\mjseqn{S(t_0) + I(t_0) + R(t_0)}}
#'   \item{`S0`}{\mjseqn{S(t_0)}}
#'   \item{`I0`}{\mjseqn{I(t_0)}}
#' }
#'
#' A warning is issued if the ODE solver cannot complete
#' the integration. A different solver may have more success
#' (e.g., consider `ode_control = list(method = "ode45")`
#' to override the default `method = "lsoda"`).
#' Using different `rtol` and `atol` might also help.
#'
#' @examples
#' # Creates a reasonable list without user input
#' par_list <- make_par_list()
#' unlist(par_list)
#'
#' # Requires time and rate parameters
#' # in units of the observation interval
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
#' @export
make_par_list <- function(dt_weeks  = 1,
                          t0        = 1000 * (365 / 7) / dt_weeks,
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
                          epsilon   = 0,
                          ode_control = list(
                            method = "lsoda",
                            rtol   = 1e-06,
                            atol   = 1e-06
                          )) {
  # Derived quantities
  gamma <- 1 / tgen
  one_year <- (365 / 7) / dt_weeks

  # If `Rnaught` was specified
  if (!is.null(Rnaught)) {
    beta_mean <- (mu / (nu * hatN0)) * Rnaught * (gamma + mu)
  # Otherwise
  } else {
    Rnaught <- ((nu * hatN0) / mu) * beta_mean / (gamma + mu)
  }


  # If `N0`, `S0`, or `I0` was not specified, then
  # obtain a value from the state of a system of
  # SIR equations after a transient of length `t0`
  if (any(is.null(c(N0, S0, I0)))) {

    # Initial state
    x_init <- c(
      S = hatN0 * (1 / Rnaught),
      logI = log(hatN0 * (1 - 1 / Rnaught) * (mu / (gamma + mu))),
      R = hatN0 * (1 - 1 / Rnaught) * (gamma / (gamma + mu))
    )

    # Time points
    t_out <- 0:t0

    # System of SIR equations
    compute_sir_rates <- function(t, y, parms) {
      y["I"] <- exp(y["logI"])
      beta <- beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
      dS <- nu * hatN0 - beta * y["S"] * y["I"] - mu * y["S"]
      dlogI <- beta * y["S"] - gamma - mu
      dR <- gamma * y["I"] - mu * y["R"]
      list(c(dS, dlogI, dR))
    }

    # List of arguments to be passed to `ode()`
    ode_args <- within(ode_control, {
      y <- x_init
      times <- t_out
      func <- compute_sir_rates
      parms <- NULL # already in environment
    })

    # Numerically integrate the system of SIR equations
    df <- as.data.frame(do.call(deSolve::ode, ode_args))

    # Assign final values of `S+I+R`, `S`, and `I`
    if (is.null(N0)) {
      N0 <- sum(df[nrow(df), c("S", "R")]) + exp(df[nrow(df), "logI"])
    }
    if (is.null(S0)) {
      S0 <- df[nrow(df), "S"]
    }
    if (is.null(I0)) {
      I0 <- exp(df[nrow(df), "logI"])
    }

    # Warn if `ode()` returned early with unrecoverable error 
    if (any(is.na(df[nrow(df), ]))) {
      warning(
        "`deSolve::ode()` could not complete the integration. ",
        "Retry with modified `ode_control`."
      )
    }
  }


  out <- as.list(environment())[names(formals(make_par_list))]
  out$ode_control <- NULL
  out
}
