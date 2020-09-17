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
#' \mjsdeqn{\out{\mathcal{R}_0 = \frac{\langle\beta\rangle N_0}{\gamma + \mu_\text{c}}}\,,}
#'
#' where \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})},
#' if `model = "sir"`, and the identity
#'
#' \mjsdeqn{\out{\mathcal{R}_0 = \frac{\langle\beta\rangle N_0}{\gamma + \mu_\text{c}} \frac{\sigma}{\sigma + \mu_\text{c}}}\,,}
#'
#' where \mjseqn{\sigma = 1 / t_\text{lat}}
#' and \mjseqn{\gamma = 1 / t_\text{inf}},
#' if `model = "seir"`.
#'
#' ## 2. Defining \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#'
#' ### SIR model
#'
#' If `model = "sir"`, then `make_par_list()` defines \mjseqn{S_0},
#' and \mjseqn{I_0} by numerically integrating the system of equations
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \mu_\text{c} N_0 - \beta(t) S I - \mu_\text{c} S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu_\text{c} I \end{align*}}
#'
#' with \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})} and
#'
#' \mjsdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{365 days}}\right) \right\rbrack}
#'
#' between times \mjseqn{t_{-(m-1)} = -(m-1) \Delta t}
#' and \mjseqn{t_0 = 0 \Delta t}, taking
#'
#' \mjsdeqn{\out{\begin{bmatrix} S(t_{-(m-1)}) \cr I(t_{-(m-1)}) \end{bmatrix} = N_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu_\text{c}}{\gamma + \mu_\text{c}} \end{bmatrix}}}
#'
#' for the initial state. This is the endemic equilibrium of the above system
#' with constant transmission rate \mjseqn{\beta(t) \equiv \langle\beta\rangle}.
#' Then `make_par_list()` assigns \mjseqn{S_0} and \mjseqn{I_0} the value of
#' \mjseqn{S(0)} and \mjseqn{I(0)}, respectively.
#'
#' ### SEIR model
#'
#' If `model = "seir"`, then `make_par_list()` defines \mjseqn{S_0},
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
#' between times \mjseqn{t_{-(m-1)} = -(m-1) \Delta t}
#' and \mjseqn{t_0 = 0 \Delta t}, taking
#'
#' \mjsdeqn{\out{\begin{bmatrix} S(t_{-(m-1)}) \cr E(t_{-(m-1)}) \cr I(t_{-(m-1)}) \end{bmatrix} = N_0 \begin{bmatrix} \frac{1}{\mathcal{R}_0} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\mu_\text{c}}{\sigma + \mu_\text{c}} \cr \big(1 - \frac{1}{\mathcal{R}_0}\big) \frac{\sigma}{\sigma + \mu_\text{c}} \frac{\mu_\text{c}}{\gamma + \mu_\text{c}} \end{bmatrix}}}
#'
#' for the initial state. This is the endemic equilibrium of the above system
#' with constant transmission rate \mjseqn{\beta(t) \equiv \langle\beta\rangle}.
#' Then `make_par_list()` assigns \mjseqn{S_0}, \mjseqn{E_0}, and \mjseqn{I_0}
#' the value of \mjseqn{S(0)}, \mjseqn{E(0)}, and \mjseqn{I(0)}, respectively.
#'
#' ## 3. Choosing \mjseqn{m}
#'
#' Choosing \mjseqn{m} such that \mjseqn{m \Delta t \sim 1000} years
#' (the default) is typically enough to ensure that \mjseqn{(S(0),I(0))}
#' (`model = "sir"`) and \mjseqn{(S(0),I(0),E(0))} (`model = "seir"`)
#' are near the attractor of the system of S(E)IR equations (specifically,
#' the state of the attractor at time \mjseqn{t_0 = 0 \Delta}). This
#' ensures that simulations by [make_data()] using these initial values
#' do not display pronounced transient dynamics, which is often desirable.
#'
#' @param dt_days \mjseqn{\lbrace\,\Delta t\,\rbrace}
#'   A numeric scalar. Observation interval in days.
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
#'   A numeric scalar. Standard deviation of the noise process
#'   realized to generate a phase shift in the seasonally forced
#'   transmission rate.
#' @param N0 \mjseqn{\lbrace\,N_0\,\rbrace}
#'   A numeric scalar. Population size at time \mjseqn{t_0 = 0 \Delta t}.
#' @param pconst \mjseqn{\lbrace\,p_\text{c}\,\rbrace}
#'   A numeric scalar. Probability that an infection is eventually reported.
#' @param model A character scalar, either `"sir"` or `"seir"`,
#'   indicating a model of disease dynamics (see Details 2).
#' @param m \mjseqn{\lbrace\,n\,\rbrace}
#'   A positive integer scalar. A system of SIR (`model = "sir"`)
#'   or SEIR (`model = "seir"`) equations is numerically
#'   integrated between times \mjseqn{t_{-(m-1)} = -(m-1) \Delta t}
#'   and \mjseqn{t_0 = 0 \Delta t} to obtain values for \mjseqn{S_0},
#'   \mjseqn{E_0}, and \mjseqn{I_0} (see Details 2 and 3).
#'
#' @return
#' A list of the arguments in the function call,
#' excluding `model` and `m` while including
#' these additional numeric scalar elements:
#'
#' \describe{
#'   \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'     Mean (long-term average) of the seasonally forced
#'     transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'   }
#'   \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'     Number of susceptible individuals at time \mjseqn{t_0 = 0 \Delta t}.
#'   }
#'   \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'     Number of exposed (infected but not infectious) individuals
#'     at time \mjseqn{t_0 = 0 \Delta t}.
#'     Included only if `model = "seir"`.
#'   }
#'   \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'     Number of infectious individuals at time \mjseqn{t_9 = 0 \Delta t}.
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
#' tinf_weeks <- 1
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
                          tlat     = 5 / dt_days,
                          tinf     = 7 / dt_days,
                          muconst  = 0.04 * dt_days / 365,
                          Rnaught  = 20,
                          alpha    = 0.08,
                          epsilon  = 0,
                          N0       = 1e06,
                          pconst   = 1,
                          model    = "sir",
                          m        = 1000 * 365 / dt_days) {
  ## Some derived quantities
  one_year <- 365 / dt_days
  if (model == "sir") {
    gamma <- 1 / (tlat + tinf)
    beta_mean <- Rnaught * (gamma + muconst) / N0
  } else if (model == "seir") {
    sigma <- 1 / tlat
    gamma <- 1 / tinf
    beta_mean <- Rnaught*(gamma+muconst)*(sigma+muconst)/(N0*sigma)
  } else {
    stop("`model` must be one of `\"sir\"` or `\"seir\"`.")
  }

  ## Time points
  m <- floor(m)
  times <- -(m-1):0

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
  S0 <- df[m, "S"]
  if (model == "seir") {
    E0 <- exp(df[m, "logE"])
  }
  I0 <- exp(df[m, "logI"])

  if (model == "sir") {
    as.list(environment())[c("dt_days", "tlat", "tinf", "muconst",
                             "Rnaught", "beta_mean", "alpha", "epsilon",
                             "N0", "S0", "I0", "pconst")]
  } else if (model == "seir") {
    as.list(environment())[c("dt_days", "tlat", "tinf", "muconst",
                             "Rnaught", "beta_mean", "alpha", "epsilon",
                             "N0", "S0", "E0", "I0", "pconst")]
  }
}
