#' \loadmathjax
#' Simulate epidemic time series data
#'
#' @description
#' `make_data()` simulates epidemic time series data using a system of
#' SIR equations and a supplied list of parameter values. Observations
#' are recorded at equally spaced time points \mjseqn{t_k = t_0 + k \Delta t}
#' (for \mjseqn{k = 0, \ldots, n}). Among other things, the simulation model
#' assumes:
#' * A seasonally forced transmission rate.
#' * Constant birth and per capita natural mortality rates.
#' * Cases reported after a fixed delay since infection,
#'   with a fixed probability.
#'
#' @param par_list A list of parameter values containing:
#'   \describe{
#'     \item{`dt_weeks`}{\mjseqn{\lbrack \Delta t \rbrack}
#'       Numeric scalar. Observation interval in weeks.
#'     }
#'     \item{`t0`}{\mjseqn{\lbrack t_0 \rbrack}
#'       Numeric scalar. Time of the first observation
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`prep`}{\mjseqn{\lbrack p_\text{rep} \rbrack}
#'       Numeric scalar. Case reporting probability.
#'     }
#'     \item{`trep`}{\mjseqn{\lbrack t_\text{rep} \rbrack}
#'       Numeric scalar. Case reporting delay in units
#'       \mjseqn{\Delta t}.
#'     }
#'     \item{`hatN0`}{\mjseqn{\lbrack \out{\widehat{N}_0} \rbrack}
#'       Numeric scalar. Population size at time \mjseqn{t = 0} years.
#'     }
#'     \item{`N0`}{\mjseqn{\lbrack N_0 \rbrack}
#'       Numeric scalar. Population size at time \mjseqn{t = t_0}.
#'     }
#'     \item{`S0`}{\mjseqn{\lbrack S_0 \rbrack}
#'       Numeric scalar. Number of susceptibles at time \mjseqn{t = t_0}.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrack I_0 \rbrack}
#'       Numeric scalar. Number of infecteds at time \mjseqn{t = t_0}.
#'     }
#'     \item{`nu`}{\mjseqn{\lbrack \nu_\text{c} \rbrack}
#'       Numeric scalar. Birth rate expressed
#'       per unit \mjseqn{\Delta t} and relative to
#'       \mjseqn{\out{\widehat{N}_0}}.
#'     }
#'     \item{`mu`}{\mjseqn{\lbrack \mu_\text{c} \rbrack}
#'       Numeric scalar. Natural mortality rate expressed
#'       per unit \mjseqn{\Delta t} and per capita.
#'     }
#'     \item{`tgen`}{\mjseqn{\lbrack t_\text{gen} \rbrack}
#'       Numeric scalar. Mean generation interval of the disease
#'       of interest in units \mjseqn{\Delta t}.
#'     }
#'     \item{`beta_mean`}{\mjseqn{\lbrack \langle\beta\rangle \rbrack}
#'       Numeric scalar. Mean (long-term average) of the seasonally
#'       forced transmission rate \mjseqn{\beta(t)} expressed per unit
#'       \mjseqn{\Delta t} per susceptible per infected.
#'     }
#'     \item{`alpha`}{\mjseqn{\lbrack \alpha \rbrack}
#'       Numeric scalar. Amplitude of the seasonally forced transmission
#'       rate \mjseqn{\beta(t)} relative to the mean.
#'     }
#'     \item{`epsilon`}{\mjseqn{\lbrack \epsilon \rbrack}
#'       Numeric scalar. Standard deviation of the random phase shift
#'       introduced to the seasonally forced transmission rate
#'       \mjseqn{\beta(t)}.
#'     }
#'   }
#' @param n Integer scalar. Time between the first and last observations
#'   in units \mjseqn{\Delta t}, so that simulated time series have `n + 1`
#'   observations. (If numeric but not integer, then `n` is replaced by
#'   `floor(n)`.)
#' @param with_dem_stoch Logical scalar. If `TRUE`, then the simulation
#'   is generated using [adaptivetau::ssa.adaptivetau()]. Otherwise, it
#'   is generated using [deSolve::ode()] (see Details).
#' @param ode_control A list of optional arguments of [deSolve::ode()],
#'   specifying options for numerical integration, such as `method`,
#'   `rtol`, and `atol`. Not used if `with_dem_stoch = TRUE` (see Details).
#'
#' @return
#' A data frame with `n + 1` rows corresponding to equally spaced times
#' \mjseqn{t_k = t_0 + k \Delta t} (for \mjseqn{k = 0, \ldots, n}), and
#' numeric columns:
#'
#' \describe{
#'   \item{`t`}{Time in units \mjseqn{\Delta t}. Equal to `par_list$t0 + 0:n`.}
#'   \item{`t_years`}{Time in years. Equal to
#'     `(par_list$t0 + 0:n) * par_list$dt_weeks * (7 / 365)`.
#'   }
#'   \item{`beta`}{Seasonally forced transmission rate expressed
#'     per unit \mjseqn{\Delta t} per susceptible per infected,
#'     without environmental noise.
#'   }
#'   \item{`beta_phi`}{Seasonally forced transmission rate expressed
#'     per unit \mjseqn{\Delta t}{dt} per susceptible per infected,
#'     with environmental noise.
#'   }
#'   \item{`N`}{Population size.}
#'   \item{`S`}{Number of susceptibles.}
#'   \item{`I`}{Number of infecteds.}
#'   \item{`Z`}{Incidence. `Z[i]` is the number of infections
#'     between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`C`}{Reported incidence. `C[i]` is the number of infections
#'     reported between times `t[i-1]` and `t[i]`, equal to the number
#'     of successes in `Z[i-round(par_list$trep)]` independent Bernoulli
#'     trials, each with success probability `par_list$prep`.
#'   }
#'   \item{`B`}{Births. `B[i]` is the number of births between times
#'     `t[i-1]` and `t[i]`.
#'   }
#'   \item{`mu`}{Per capita natural mortality rate. `mu[i]` is rate
#'     at time `t[i]`. Equal to `rep(par_list$mu, n + 1)` and included
#'     for convenience, mainly because [estimate_beta_si()] requires
#'     a data frame with columns `Z`, `B`, and `mu`.
#'   }
#' }
#'
#' The data frame has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(make_data, arg_list)`,
#' preceded by a call to `set.seed()`.
#'
#' @details
#' # Details
#' 
#' ## Simulation model
#' 
#' `make_data()` simulates epidemic time series data using the system of
#' SIR equations below, which includes a fourth equation for cumulative
#' births and a fifth equation for cumulative incidence:
#'
#' \mjseqn{\out{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu \widehat{N}_0 - \beta_\phi(t) S I - \mu S\,, \cr \frac{\text{d}I}{\text{d}t} &= \beta_\phi(t) S I - \gamma I - \mu I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu R\,, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu \widehat{N}_0\,, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta_\phi(t) S I\,. \end{align}}}
#' 
#' Here, \mjseqn{\gamma = 1 / t_\text{gen}},
#' 
#' \mjsdeqn{\beta_\phi(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}} + \phi(t;\epsilon)\right) \right\rbrack\,,}
#' 
#' and \mjseqn{\phi(t;\epsilon)} is the linear interpolant of noise
#' \mjseqn{\lbrace(t_k;\Phi_k)\rbrace_{k=0}^n} with
#' 
#' \mjsdeqn{\Phi_k \sim \mathrm{Normal}(0,\epsilon^2)\,,}
#'
#' modeling **environmental stochasticity**.
#'
#' `make_data()` generates observations of the above system at equally spaced
#' times \mjseqn{t_k = t_0 + k \Delta t} (for \mjseqn{k = 0, \ldots, n}) by
#' either
#' (i) numerically integrating the ODE using [deSolve::ode()]
#' (`with_dem_stoch = FALSE`), or
#' (ii) realizing a corresponding continuous-time stochastic process using
#' [adaptivetau::ssa.adaptivetau()] (`with_dem_stoch = TRUE`).
#' Both methods use initial state
#'
#' \mjsdeqn{\begin{bmatrix} S(t_0) \cr I(t_0) \cr R(t_0) \cr B_\text{cum}(t_0) \cr Z_\text{cum}(t_0) \end{bmatrix} = \begin{bmatrix} S_0 \cr I_0 \cr N_0 - S_0 - I_0 \cr 0 \cr 0 \end{bmatrix}\,.}
#' 
#' The latter method defines event probabilities as proportional to
#' terms in the ODE and models **demographic stochasticity**. Disease
#' fadeout in simulations with demographic stochasticity is prevented
#' by setting the probability of infected decrease to zero whenever
#' the number of infecteds is one.
#'
#' `make_data()` calculates births \mjseqn{B} and incidence \mjseqn{Z}
#' from cumulative births \mjseqn{B_\text{cum}} and cumulative incidence
#' \mjseqn{Z_\text{cum}} via first differences,
#'
#' \mjsdeqn{\begin{align} B(t_k) &= B_\text{cum}(t_k) - B_\text{cum}(t_{k-1}), \cr Z(t_k) &= Z_\text{cum}(t_k) - Z_\text{cum}(t_{k-1}), \end{align}}
#'
#' then simulates reported incidence \mjseqn{C} from
#' incidence \mjseqn{Z} via lagged binomial sampling,
#'
#' \mjsdeqn{C(t_{k+r}) \sim \mathrm{Binomial}\big(Z(t_k),p_\text{rep}\big)\,,}
#'
#' modeling **observation error** (random under-reporting of cases with
#' a delay between infection and reporting). Here,
#' \mjseqn{r = \mathrm{nint}(t_\text{rep} / \Delta t)}.
#'
#' @examples
#' # Deterministic simulation
#' par_list <- make_par_list(dt_weeks = 1)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7, # number of weeks in 20 years
#'   with_dem_stoch = FALSE
#' )
#' head(df)
#' 
#' # Stochastic simulation
#' par_list <- make_par_list(
#'   dt_weeks = 1,
#'   epsilon = 0.5, # environmental stochasticity
#'   prep = 0.5 # random under-reporting of cases
#' )
#' set.seed(2146)
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7,
#'   with_dem_stoch = TRUE # demographic stochasticity
#' )
#' head(df)
#'
#' @export
make_data <- function(par_list       = make_par_list(),
                      n              = 20 * 365 / 7,
                      with_dem_stoch = TRUE,
                      ode_control    = list(
                        method = "lsoda",
                        rtol   = 1e-06,
                        atol   = 1e-06
                      )) {
  ## 1. Set-up -----------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(
    par_list[c("dt_weeks", "t0", "prep", "trep", "hatN0", "N0", "S0", "I0",
               "nu", "mu", "tgen", "beta_mean", "alpha", "epsilon")],
    envir = environment()
  )

  # Derived quantities
  gamma <- 1 / tgen
  one_year <- (365 / 7) / dt_weeks

  # Time points
  n <- floor(n)
  t_out <- t0 + 0:n

  # Environmental noise
  phi <- stats::rnorm(
    n    = length(t_out),
    mean = 0,
    sd   = epsilon
  )

  # Linear interpolant of noise
  interpolate_phi <- stats::approxfun(
    x      = t_out,
    y      = phi,
    method = "linear",
    rule   = 2 # return `y[1]` and `y[length(y)]` outside range of `x`
  )

  # Seasonally forced transmission rate
  # without environmental noise
  beta <- function(t) {
    beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
  }
  # with environmental noise
  beta_phi <- function(t) {
    phi <- interpolate_phi(t)
    beta_mean * (1 + alpha * cos(2 * pi * t / one_year + phi))
  }


  ## 2.(a) Simulate SIR equations ... ------------------------------------
  ##       if with demographic stochasticity

  if (with_dem_stoch) {

    ## NOTE: The adaptivetau package insists that simulations start at
    ##       time 0. To get simulations from time `t0`, we must take care
    ##       to add `t0` to the simulation time `t` where necessary.

    # Initial state
    x_init <- c(
      S = ceiling(S0),                             # susceptibles
      I = ceiling(I0),                             # infecteds
      R = ceiling(N0) - ceiling(S0) - ceiling(I0), # removeds
      Bcum = 0,                                    # cumulative births
      Zcum = 0                                     # cumulative incidence
    )

    # Transition events
    event_list <- list(
      c(S = 1, Bcum = 1),         # birth
      c(S = -1, I = 1, Zcum = 1), # infection
      c(I = -1, R = 1),           # removal
      c(S = -1),                  # natural mortality
      c(I = -1),
      c(R = -1)
    )

    # Transition event rates
    compute_event_rates <- function(x, params, t) {
      c(
        nu * hatN0,                         # birth
        beta_phi(t + t0) * x["S"] * x["I"], # infection
        gamma * x["I"] * (x["I"] > 1),      # removal
        mu * x["S"],                        # natural mortality
        mu * x["I"] * (x["I"] > 1),
        mu * x["R"]
      )
    }

    # Generate a realization of the stochastic process
    df <- as.data.frame(adaptivetau::ssa.adaptivetau(
      x_init, event_list, compute_event_rates,
      params    = NULL, # already in environment
      tf        = n,    # final time point
      tl.params = list( # other instructions:
        epsilon     = 0.05,
        delta       = 0.05,
        maxtau      = 0.5, # adaptive time step must not exceed 1
        extraChecks = TRUE
      )
    ))
    colnames(df) <- c("t", "S", "I", "R", "Bcum", "Zcum")
    df$t <- df$t + t0

    ## NOTE: `ssa.adaptivetau()` returns time series with unequal spacing
    ##       (time step varies between 0 and `maxtau`), but we desire
    ##       equal spacing (time step equal to 1, as in `t_out`)

    # For each time in `t_out`, find the index
    # of the last previous time in `df$t`
    ind_state_out <- sapply(t_out, function(t) max(which(df$t <= t)))

    # Take the state at those times to be the state at times `t_out`
    df <- df[ind_state_out, ]
    df$t <- t_out
    rownames(df) <- NULL


  ## 2.(b) Simulate SIR equations ... ------------------------------------
  ##       if without demographic stochasticity

  } else {

    # Initial state
    x_init <- c(
      S    = S0,           # susceptibles
      logI = log(I0),      # log infecteds
      R    = N0 - S0 - I0, # removeds
      Bcum = 0,            # cumulative births
      Zcum = 0             # cumulative incidence
    )

    # System of SIR equations with additional equations for
    # cumulative births `Bcum` and cumulative incidence `Zcum`
    compute_sir_rates <- function(t, y, parms) {
      y["I"] <- exp(y["logI"])
      beta_phi <- beta_phi(t)
      dBcum <- nu * hatN0
      dZcum <- beta_phi * y["S"] * y["I"]
      dS <- dBcum - dZcum - mu * y["S"]
      dlogI <- beta_phi * y["S"] - gamma - mu
      dR <- gamma * y["I"] - mu * y["R"]
      list(c(dS, dlogI, dR, dBcum, dZcum))
    }

    # Create a list of arguments to be passed to `ode()`
    ode_args <- within(ode_control, {
      y <- x_init
      times <- t_out
      func <- compute_sir_rates
      parms <- NULL # already in environment
    })

    # Numerically integrate the system of SIR equations
    df <- as.data.frame(do.call(deSolve::ode, ode_args))
    colnames(df) <- c("t", "S", "logI", "R", "Bcum", "Zcum")
    df$I <- exp(df$logI)

    # Warn if `ode()` returned early with unrecoverable error.
    # Append rows of `NA` until `nrow(df) = length(t_out)`.
    if (any(is.na(df[nrow(df), ]))) {
      warning(
        "`deSolve::ode()` could not complete the integration. ",
        "Retry with modified `ode_control`."
      )
    }  
    if (nrow(df) < length(t_out)) {
      df[(nrow(df)+1):length(t_out), ] <- NA
      df$t <- t_out
    }

  }


  ## 3. Append other desired variables -----------------------------------

  df <- transform(df,
    t_years  = t * dt_weeks * (7 / 365),
    beta     = beta(t),
    beta_phi = beta_phi(t),
    N        = S + I + R,
    Z        = c(NA, diff(Zcum)), # `Z` from `Zcum` via first differences
    B        = c(NA, diff(Bcum)), # `B` from `Bcum` via first differences
    mu       = mu
  )

  # `C` from `Z` via lagged binomial sampling. `rbinom()`
  # warning about `NA` in `size` argument is safe to suppress.
  df$C <- suppressWarnings({
    stats::rbinom(
      n    = nrow(df),    # number of experiments
      size = round(df$Z), # number of Bernoulli trials
      p    = prep         # success probability
    )
  })
  trepr <- round(trep)
  df$C <- c(rep(NA, trepr), df$C[1:(nrow(df)-trepr)])


  df <- df[, c("t", "t_years", "beta", "beta_phi",
               "N", "S", "I", "Z", "C", "B", "mu")]
  attr(df, "call") <- match.call()
  attr(df, "arg_list") <- arg_list
  df
}

# For `R CMD check`
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("dt_weeks", "t0", "prep", "trep", "hatN0", "N0", "S0", "I0",
      "nu", "mu", "tgen", "beta_mean", "alpha", "epsilon", "S", "R",
      "Zcum", "Bcum")
  )
}
