#' \loadmathjax
#' Create a list of parameter values for simulations
#'
#' @description
#' Creates a list of parameter values that can be assigned
#' to the argument `par_list` of [make_data()] if defaults
#' are accepted for its arguments `mu`, `nu`, `beta`, and
#' `p`.
#'
#' @details
#' # Details
#'
#' ## 1. Defining \mjseqn{\langle\beta\rangle}
#'
#' `make_par_list()` defines \mjseqn{\langle\beta\rangle} as a function
#' of \mjseqn{\mathcal{R}_0} by enforcing the identity
#'
#' \mjsdeqn{\mathcal{R}_0 = \frac{\langle\beta\rangle N_0}{\gamma + \mu}}
#'
#' if `model = "sir"`, where \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})},
#' and the identity
#'
#' \mjsdeqn{\mathcal{R}_0 = \frac{\langle\beta\rangle N_0}{\gamma + \mu} \frac{\sigma}{\sigma + \mu}}
#'
#' if `model = "seir"`, where \mjseqn{\sigma = 1 / t_\text{lat}}
#' and \mjseqn{\gamma = 1 / t_\text{inf}}.
#'
#' ## 2. Defining \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#'
#' ### SIR model
#'
#' If `model == "sir"`, then `make_par_list()` defines \mjseqn{S_0},
#' and \mjseqn{I_0} by numerically integrating the system of equations
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \mu_\text{c} N_0 - \beta(t) S I - \mu_\text{c} S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu_\text{c} I \end{align*}}
#'
#' with \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})} and
#'
#' \mjsdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{365 days}}\right) \right\rbrack}
#'
#' between times \mjseqn{t = t_{-(n-1)} = -(n-1) \Delta t}
#' and \mjseqn{t = 0} years, taking
#'
#' \mjsdeqn{\begin{bmatrix} S(t_{-(n-1)}) \cr I(t_{-(n-1)}) \end{bmatrix} = N_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu}{\gamma + \mu} \end{bmatrix}}
#'
#' for the initial state. This is the endemic equilibrium of the above system
#' with constant transmission rate \mjseqn{\beta(t) \equiv \langle\beta\rangle}.
#' Then `make_par_list()` assigns \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#' the value of \mjseqn{S(0)}, \mjseqn{E(0)}, and \mjseqn{I(0)}, respectively.
#'
#' ### SEIR model
#'
#' If `model == "seir"`, then `make_par_list()` defines \mjseqn{S_0},
#' \mjseqn{E_0}, and \mjseqn{I_0} by numerically integrating the system
#' of equations
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \mu_\text{c} N_0 - \beta(t) S I - \mu_\text{c} S \cr \frac{\text{d}E}{\text{d}t} &= \beta(t) S I - \sigma E - \mu_\text{c} E \cr \frac{\text{d}I}{\text{d}t} &= \sigma E - \gamma I - \mu_\text{c} I \end{align*}}
#'
#' with \mjseqn{\sigma = 1 / t_\text{lat}}, \mjseqn{\gamma = 1 / t_\text{inf}},
#' and
#'
#' \mjsdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{365 days}}\right) \right\rbrack}
#'
#' between times \mjseqn{t = t_{-(n-1)} = -(n-1) \Delta t}
#' and \mjseqn{t = 0} years, taking
#'
#' \mjsdeqn{\begin{bmatrix} S(t_{-(n-1)}) \cr E(t_{-(n-1)}) \cr I(t_{-(n-1)}) \end{bmatrix} = N_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu}{\sigma + \mu} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\sigma}{\sigma + \mu} \frac{\mu}{\gamma + \mu} \end{bmatrix}}
#'
#' for the initial state. This is the endemic equilibrium of the above system
#' with constant transmission rate \mjseqn{\beta(t) \equiv \langle\beta\rangle}.
#' Then `make_par_list()` assigns \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#' the value of \mjseqn{S(0)}, \mjseqn{E(0)}, and \mjseqn{I(0)}, respectively.
#'
#' ## 3. Choosing \mjseqn{n}
#'
#' Choosing \mjseqn{n} such that \mjseqn{n \Delta t \sim 1000}
#' years (the default) is typically enough to ensure that
#' \mjseqn{\big(S(0),I(0)\big)} (`model = "sir"`) or
#' \mjseqn{\big(S(0),I(0),E(0)\big)} (`model = "seir"`)
#' are near the attractor of the system of S(E)IR equations,
#' which can be desirable when using [make_data()] to simulate
#' epidemic time series.
#'
#' @param dt_days \mjseqn{\lbrace\,\Delta t\,\rbrace}
#'   A numeric scalar. Observation interval in days.
#' @param N0 \mjseqn{\lbrace\,N_0\,\rbrace}
#'   A numeric scalar. Population size at time \mjseqn{t = 0} years.
#' @param tlat \mjseqn{\lbrace\,t_\text{lat}/\Delta t\,\rbrace}
#'   A numeric scalar. Mean latent period
#'   of the disease of interest in units \mjseqn{\Delta t}.
#' @param tinf \mjseqn{\lbrace\,t_\text{inf}/\Delta t\,\rbrace}
#'   A numeric scalar. Mean infectious period
#'   of the disease of interest in units \mjseqn{\Delta t}.
#' @param muconst \mjseqn{\lbrace\,\mu_\text{c} \Delta t\,\rbrace}
#'   A numeric scalar. Per capita natural mortality rate expressed
#'   per unit \mjseqn{\Delta t}.
#' @param Rnaught \mjseqn{\lbrace\,\mathcal{R}_0\,\rbrace}
#'   A numeric scalar. Basic reproduction number
#'   of the disease of interest.
#' @param alpha \mjseqn{\lbrace\,\alpha\,\rbrace}
#'   A numeric scalar. Amplitude of the seasonally forced
#'   transmission rate relative to its mean.
#' @param epsilon \mjseqn{\lbrace\,\epsilon\,\rbrace}
#'   A numeric scalar. Standard deviation of the standard
#'   normally distributed phase shift in the seasonally forced
#'   transmission rate.
#' @param pconst \mjseqn{\lbrace\,p_\text{c}\,\rbrace}
#'   A numeric scalar. Probability that an infection is eventually reported.
#' @param model A character scalar, either `"sir"` or `"seir"`,
#'   indicating a system of equations to be numerically integrated
#'   (see Details 2).
#' @param n \mjseqn{\lbrace\,n\,\rbrace}
#'   A positive integer scalar. A system of SIR (`model = "sir"`)
#'   or SEIR (`model = "seir"`) equations will be numerically integrated
#'   between times \mjseqn{t = -(n-1) \Delta t} and \mjseqn{t = 0} years
#'   to obtain values for \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#'   (see Details 2 and 3).
#'
#' @return
#' A list of the arguments in the function call,
#' excluding `model` and `n` while including
#' these additional numeric scalar elements:
#'
#' \describe{
#'   \item{`tgen`}{\mjseqn{\lbrace\,t_\text{gen} / \Delta t\,\rbrace}
#'     Mean generation interval of the disease of interest
#'     in units \mjseqn{\Delta t}. Equal to `tlat + tinf`.
#'   }
#'   \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'     Mean (long-term average) of the seasonally forced
#'     transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'   }
#'   \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'     Number of susceptible individuals at time \mjseqn{t = 0} years.
#'   }
#'   \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'     Number of exposed (infected but not infectious) individuals
#'     at time \mjseqn{t = 0} years. Included only if `model = "seir"`.
#'   }
#'   \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'     Number of infectious individuals at time \mjseqn{t = 0} years.
#'   }
#' }
#'
#' This list can be assigned to the argument `par_list` of `[make_data()]`
#' if defaults are accepted for its arguments `mu`, `nu`, `beta`, and `p`.
#' In this case, care must be taken to ensure that the argument `model`
#' of [make_data()] matches that of `make_par_list()`.
#'
#' @examples
#' ## Creates a plausible list for measles by default
#' pl <- make_par_list()
#' unlist(pl)
#'
#' ## Time (rate) parameters should be specified
#' ## in units (per unit) of the observation interval
#' dt_days <- 7
#' tlat_days <- 5
#' tinf_weeks <- 8 / 7
#' muconst_peryear <- 0.04
#' pl <- make_par_list(
#'   dt_days = dt_days,
#'   tlat    = tlat_days / dt_days,
#'   tinf    = tinf_weeks * 7 / dt_days,
#'   muconst = muconst_peryear * dt_days / 365
#' )
#' unlist(pl)
#'
#' @seealso [make_data()]
#' @export
#' @importFrom deSolve ode
make_par_list <- function(dt_days  = 7,
                          N0       = 1e06,
                          tlat     = 5 / dt_days,
                          tinf     = 7 / dt_days,
                          muconst  = 0.04 * dt_days / 365,
                          Rnaught  = 20,
                          alpha    = 0.08,
                          epsilon  = 0,
                          pconst   = 1,
                          model    = "sir",
                          n        = 1000 * 365 / dt_days) {
  ## Some derived quantities
  one_year <- 365 / dt_days
  tgen <- tlat + tinf
  if (model == "sir") {
    gamma <- 1 / tgen
    beta_mean <- Rnaught * (gamma + muconst) / N0
  } else if (model == "seir") {
    sigma <- 1 / tlat
    gamma <- 1 / tinf
    beta_mean <- Rnaught*(gamma+muconst)*(sigma+muconst)/(N0*sigma)
  } else {
    stop("`model` must be one of `\"sir\"` or `\"seir\"`.")
  }

  ## Time points
  n <- floor(n)
  times <- -(n - 1):0

  if (model == "sir") {

    ## Initial state
    x_init <- c(
      S    = N0 * (1 / Rnaught),
      logI = log(N0 * (1 - 1 / Rnaught) * (muconst / (gamma + muconst)))
    )

    ## System of SIR equations
    compute_ode_rates <- function(t, y, parms) {
      y["I"] <- exp(y["logI"])
      beta <- beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
      dS <- muconst * N0 - beta * y["S"] * y["I"] - muconst * y["S"]
      dlogI <- beta * y["S"] - gamma - muconst
      list(c(dS, dlogI))
    }

  } else if (model == "seir") {

    ## Initial state
    x_init <- c(
      S    = N0*(1/Rnaught),
      logE = log(N0*(1-1/Rnaught)*(muconst/(sigma+muconst))),
      logI = log(N0*(1-1/Rnaught)*(sigma/(sigma+muconst))*(muconst/(gamma+muconst)))
    )

    ## System of SEIR equations
    compute_ode_rates <- function(t, y, parms) {
      y["E"] <- exp(y["logE"])
      y["I"] <- exp(y["logI"])
      beta <- beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
      dS <- muconst * N0 - beta * y["S"] * y["I"] - muconst * y["S"]
      dlogE <- beta * y["S"] * y["I"] / y["E"] - sigma - muconst
      dlogI <- sigma * y["E"] / y["I"] - gamma - muconst
      list(c(dS, dlogE, dlogI))
    }

  }

  ## Numerically integrate the system of S(E)IR equations
  df <- as.data.frame(
    ode(
      y     = x_init,
      times = times,
      func  = compute_ode_rates,
      parms = NULL, # found in enclosing environment of `compute_ode_rates()`
      hmax  = 1     # avoids error when `length(times) = 1`
    )
  )

  ## Assign final values of `S`, `E`, and `I`
  S0 <- df[n, "S"]
  if (model == "seir") {
    E0 <- exp(df[n, "logE"])
  }
  I0 <- exp(df[n, "logI"])

  if (model == "sir") {
    as.list(environment())[c("dt_days", "N0", "S0", "I0",
                             "tlat", "tinf", "tgen", "muconst", "Rnaught",
                             "beta_mean", "alpha", "epsilon", "pconst")]
  } else if (model == "seir") {
    as.list(environment())[c("dt_days", "N0", "S0", "E0", "I0",
                             "tlat", "tinf", "tgen", "muconst", "Rnaught",
                             "beta_mean", "alpha", "epsilon", "pconst")]
  }
}
