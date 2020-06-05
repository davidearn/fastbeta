#' \loadmathjax
#' Create a list of parameter values for simulations
#'
#' @description
#' `make_par_list()` creates a list of parameter values to be passed
#' as an argument of [make_data()], which simulates epidemic time series
#' data using an SIR model.
#'
#' @section Concordance of \mjeqn{\mathcal{R}_0}{calR_0} and \mjeqn{\langle\beta\rangle}{<beta>}:
#' `make_par_list()` enforces the identity
#'
#' \mjdeqn{\mathcal{R}_0 = \frac{\nu \widehat{N}_0}{\mu} \cdot \frac{\langle\beta\rangle}{\gamma + \mu},}{calR_0 = ((nu hatN_0) / mu) (<beta> / (gamma + mu)),}
#'
#' where \mjeqn{\gamma = 1 / t_\text{gen}}{gamma = 1 / t_gen}.
#' This means that it calculates \mjeqn{\langle\beta\rangle}{<beta>}
#' as a function of \mjeqn{\mathcal{R}_0}{calR_0}, unless the user
#' explicitly sets `Rnaught = NULL`, in which case it calculates
#' \mjeqn{\mathcal{R}_0}{calR_0} as a function of
#' \mjeqn{\langle\beta\rangle}{<beta>}. Hence argument `beta_mean`
#' is optional if `Rnaught` is specified, but mandatory otherwise.
#' 
#' @section Missing \mjseqn{N_0}, \mjseqn{S_0}, or \mjseqn{I_0}:
#' If \mjseqn{N_0}, \mjseqn{S_0}, or \mjseqn{I_0} (`N0`, `S0`, or `I0`)
#' is not specified in the function call, then, with a call to
#' [deSolve::ode()], `make_par_list()` numerically integrates the
#' system of SIR equations
#'
#' \mjdeqn{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu \widehat{N}_0 - \beta(t) S I - \mu S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu I \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu R \end{align}}{(1) dS/dt = nu hatN_0 - beta(t) S I - mu S    (2) dI/dt = beta(t) S I - gamma I - mu I    (3) dR/dt = gamma I - mu R}
#' 
#' with \mjeqn{\gamma = 1 / t_\text{gen}}{gamma = 1 / t_gen} and
#' 
#' \mjdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}}\right) \right\rbrack}{beta(t) = <beta> (1 + alpha cos(2 pi t / (1 year)))}
#' 
#' between times \mjseqn{t = 0} years and \mjseqn{t = t_0}, taking
#'
#' \mjdeqn{\begin{bmatrix} S(0) \cr I(0) \cr R(0) \end{bmatrix} = \widehat{N}_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu}{\gamma + \mu} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\gamma}{\gamma + \mu} \end{bmatrix},}{(S(0),I(0),R(0)) = (hatS,hatI,hatR),}
#'
#' for the initial state. This is the endemic equilibrium
#' of the above system with constant transmission rate
#' \mjeqn{\beta(t) \equiv \langle\beta\rangle}{beta(t) = <beta>},
#' scaled by \mjeqn{\mu/\nu}{mu/nu}
#' to enforce
#'
#' \mjdeqn{S(0) + I(0) + R(0) = \widehat{N}_0}{S(0) + I(0) + R(0) = hatN_0}
#'
#' as the initial population size. Then `make_par_list()` assigns `N0`,
#' `S0`, and `I0` (only those not specified in the function call) the
#' values
#'
#' \describe{
#'   \item{`N0`}{\mjseqn{S(t_0) + I(t_0) + R(t_0).}}
#'   \item{`S0`}{\mjseqn{S(t_0).}}
#'   \item{`I0`}{\mjseqn{I(t_0).}}
#' }
#'
#' A warning is issued if the ODE solver cannot complete
#' the integration. A different solver may have more success
#' (e.g., consider `ode_control = list(method = "ode45")`
#' to override the default `method = "lsoda"`).
#' Using different `rtol` and `atol` might also help.
#'
#' @param dt_weeks \mjeqn{\lbrack \Delta t \rbrack}{\[dt\]}
#'   Numeric scalar. Observation interval in weeks.
#' @param t0 \mjeqn{\lbrack t_0 \rbrack}{\[t_0\]}
#'   Numeric scalar. Time of the first observation
#'   in units \mjeqn{\Delta t}{dt}.
#' @param prep \mjeqn{\lbrack p_\text{rep} \rbrack}{\[p_rep\]}
#'   Numeric scalar. Case reporting probability.
#' @param trep \mjeqn{\lbrack t_\text{rep} \rbrack}{\[t_rep\]}
#'   Numeric scalar. Case reporting delay in units \mjeqn{\Delta t}{dt}.
#' @param hatN0 \mjeqn{\lbrack \widehat{N}_0 \rbrack}{\[hatN_0\]}
#'   Numeric scalar. Population size at time \mjseqn{t = 0} years
#'   (see Details).
#' @param N0 \mjeqn{\lbrack N_0 \rbrack}{\[N_0\]}
#'   Numeric scalar. Population size at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details).
#' @param S0 \mjeqn{\lbrack S_0 \rbrack}{\[S_0\]}
#'   Numeric scalar. Number of susceptibles at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details).
#' @param I0 \mjeqn{\lbrack I_0 \rbrack}{\[I_0\]}
#'   Numeric scalar. Number of infecteds at time \mjseqn{t = t_0}.
#'   This argument is optional (see Details).
#' @param nu \mjeqn{\lbrack \nu \rbrack}{\[nu\]}
#'   Numeric scalar. Birth rate expressed per unit \mjeqn{\Delta t}{dt}
#'   and relative to \mjeqn{\widehat{N}_0}{hatN_0}.
#' @param mu \mjeqn{\lbrack \mu \rbrack}{\[mu\]}
#'   Numeric scalar. Natural mortality rate expressed per unit
#'   \mjeqn{\Delta t}{dt} and per capita.
#' @param tgen \mjeqn{\lbrack t_\text{gen} \rbrack}{\[t_gen\]}
#'   Numeric scalar. Mean generation interval of the disease
#'   of interest in units \mjeqn{\Delta t}{dt}.
#' @param Rnaught \mjeqn{\lbrack \mathcal{R}_0 \rbrack}{\[calR_0\]}
#'   Numeric scalar. Basic reproduction number of the disease
#'   of interest. If specified (not `NULL`), then `beta_mean` is
#'   calculated internally as a function of `Rnaught` (see Details).
#' @param beta_mean \mjeqn{\lbrack \langle\beta\rangle \rbrack}{\[<beta>\]}
#'   Numeric scalar. Mean (long-term average) of the seasonally
#'   forced transmission rate \mjeqn{\beta(t)}{beta(t)} expressed
#'   per unit \mjeqn{\Delta t}{dt} per susceptible per infected.
#'   Ignored if `Rnaught` is specified (not `NULL`) (see Details).
#' @param alpha \mjeqn{\lbrack \alpha \rbrack}{\[alpha\]}
#'   Numeric scalar. Amplitude of the seasonally forced transmission
#'   rate \mjeqn{\beta(t)}{beta(t)} relative to the mean.
#' @param epsilon \mjeqn{\lbrack \epsilon \rbrack}{\[epsilon\]}
#'   Numeric scalar. Standard deviation of the random phase shift
#'   introduced to the seasonally forced transmission rate
#'   \mjeqn{\beta(t)}{beta(t)}.
#' @param ode_control A list of optional arguments of [deSolve::ode()],
#'   specifying options for numerical integration, such as `method`,
#'   `rtol`, and `atol` (see Details).
#'
#' @return
#' A list of the arguments of `make_par_list()` (excluding
#' `ode_control`) with values for `Rnaught`, `beta_mean`,
#' `N0`, `S0`, and `I0` if not supplied in the function call
#' (see Details).
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
if (any(is.na(c(N0, S0, I0)))) {

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
  mat <- do.call(deSolve::ode, ode_args)

  # Assign final values of `S+I+R`, `S`, and `I`
  if (is.null(N0)) {
    N0 <- sum(mat[nrow(mat), c("S", "R")]) + exp(mat[nrow(mat), "logI"])
  }
  if (is.null(S0)) {
    S0 <- mat[nrow(mat), "S"]
  }
  if (is.null(I0)) {
    I0 <- exp(mat[nrow(mat), "logI"])
  }

  # Warn if `ode()` returned early with unrecoverable error 
  if (any(is.na(mat[nrow(mat), ]))) {
    warning(
      "`deSolve::ode()` could not complete the integration. ",
      "Retry with modified `ode_control`.",
      call. = FALSE
    )
  }
}


out <- as.list(environment())[names(formals(make_par_list))]
out$ode_control <- NULL
out
}
