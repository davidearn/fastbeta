#' \loadmathjax
#' Simulate epidemic time series data
#'
#' @description
#' Simulates epidemic time series data using a system of S(E)IR equations
#' (see Details 1). Observations are recorded at equally spaced time points
#' \mjseqn{t_i = t_0 + i \Delta t}. Users can specify any time-varying birth,
#' death, and transmission rates and any discrete distribution of the time
#' from infection to reporting. By default, the vital rates are constant and
#' the transmission rate is seasonally forced (see Details 2).
#'
#' @details
#' # Details
#'
#' ## 1. Simulation model
#'
#' `make_data()` simulates epidemic time series data using a system
#' of SIR or SEIR equations, which includes additional equations for
#' cumulative births and cumulative incidence. If `model = "sir"`,
#' then the system is
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) N_0 - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu(t) N_0\,, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta(t) S I\,, \end{align*}}
#'
#' where \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})}.
#' If `model = "seir"`, then the system is
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) N_0 - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}E}{\text{d}t} &= \beta(t) S I - \sigma E - \mu(t) E\,, \cr \frac{\text{d}I}{\text{d}t} &= \sigma E - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu(t) N_0\,, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta(t) S I\,, \end{align*}}
#'
#' where \mjseqn{\sigma = 1 / t_\text{lat}} and
#' \mjseqn{\gamma = 1 / t_\text{inf}}. Note that setting
#' \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})}
#' in the SIR model guarantees that the implied mean
#' generation interval matches that of the corresponding
#' SEIR model.
#'
#' `make_data()` generates observations of the system at equally spaced
#' times \mjseqn{t_i = i \Delta t} (for \mjseqn{i = 0, \ldots, n-1}) by
#' either
#' (i) numerically integrating the ODE using [deSolve::ode()]
#' (`with_ds = FALSE`), or
#' (ii) realizing a corresponding continuous-time stochastic process
#' using [adaptivetau::ssa.adaptivetau()] (`with_ds = TRUE`).
#' Both methods (`with_ds = FALSE` and `with_ds = TRUE`) use
#' the initial state
#'
#' \mjsdeqn{\begin{bmatrix} S(0) \cr I(0) \cr R(0) \cr B_\text{cum}(0) \cr Z_\text{cum}(0) \end{bmatrix} = \begin{bmatrix} S_0 \cr I_0 \cr N_0 - S_0 - I_0 \cr 0 \cr 0 \end{bmatrix}}
#'
#' if `model = "sir"` or
#'
#' \mjsdeqn{\begin{bmatrix} S(0) \cr E(0) \cr I(0) \cr R(0) \cr B_\text{cum}(0) \cr Z_\text{cum}(0) \end{bmatrix} = \begin{bmatrix} S_0 \cr E_0 \cr I_0 \cr N_0 - S_0 - E_0 - I_0 \cr 0 \cr 0 \end{bmatrix}}
#'
#' if `model = "seir"`.
#'
#' The latter method (`with_ds = TRUE`) defines event probabilities as
#' proportional to terms in the ODE, modeling **demographic stochasticity**.
#' Disease fadeout in simulations with demographic stochasticity is prevented
#' by setting the probability of a decrease in the infected population size
#' to zero whenever the infected population size is one.
#'
#' `make_data()` calculates births \mjseqn{B} and incidence \mjseqn{Z}
#' from cumulative births \mjseqn{B_\text{cum}} and cumulative incidence
#' \mjseqn{Z_\text{cum}} via first differences:
#'
#' \mjsdeqn{\begin{align*} B(t_i) &= B_\text{cum}(t_i) - B_\text{cum}(t_{i-1}), \cr Z(t_i) &= Z_\text{cum}(t_i) - Z_\text{cum}(t_{i-1}). \end{align*}}
#'
#' It then does binomial sampling to generate, for each observation
#' of \mjseqn{Z}, a number of infections that is eventually reported:
#'
#' \mjsdeqn{Z_\text{rep}(t_i) \sim \mathrm{Binomial}\big(Z(t_i),p_\text{rep}\big)\,.}
#'
#' Finally, it generates reported incidence \mjseqn{C} by sampling from
#' a supplied reporting delay distribution and binning infections counted
#' by \mjseqn{Z_\text{rep}} forward in time. Compared to \mjseqn{Z},
#' \mjseqn{C} models **observation error**, i.e., under-reporting of
#' cases with a delay between infection and reporting.
#'
#' ## 2. Default parametrization
#'
#' By default, \mjseqn{\mu} and and \mjseqn{\nu} are modeled as constant,
#' while \mjseqn{\beta} is modeled by a seasonal forcing function with
#' environmental noise. Specifically,
#'
#' \mjsdeqn{\begin{align*} \mu(t) &= \mu_\text{c}\,, \cr \nu(t) &= \mu_\text{c}\,, \cr \beta(t) &= \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{365 days}} + \phi(t;\epsilon^2)\right) \right\rbrack\,, \end{align*}}
#'
#' where \mjseqn{\phi(t;\epsilon^2)} is the linear interpolant of noise
#' \mjseqn{\lbrace(t_i;\Phi_i)\rbrace_{i=0}^{n-1}} with
#'
#' \mjsdeqn{\Phi_i \sim \mathrm{Normal}(0,\epsilon^2)\,.}
#'
#' Care is taken to ensure that random number generation occurs in the
#' enclosing environment rather than the body of the function assigned
#' by default to `beta`. This ensures that the function has noisy but
#' non-random output, which is helpful for reproducibility. For details,
#' see [make_beta()].
#'
#' The default specification of `mu`, `nu`, and `beta` requires
#' that `par_list` contains these additional numeric scalar elements:
#'
#'   \describe{
#'     \item{`muconst`}{\mjseqn{\lbrace\,\mu_\text{c} \Delta t\,\rbrace}
#'       Birth rate (relative to \mjseqn{N_0}) expressed per unit
#'       \mjseqn{\Delta t} *and* per capita natural mortality rate
#'       expressed per unit \mjseqn{\Delta t}.
#'     }
#'     \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'       Mean (long-term average) of the seasonally forced
#'       transmission rate expressed per unit \mjseqn{\Delta t}
#'       per susceptible individual per infectious individual.
#'     }
#'     \item{`alpha`}{\mjseqn{\lbrace\,\alpha\,\rbrace}
#'       Amplitude of the seasonally forced transmission rate
#'       relative to the mean.
#'     }
#'     \item{`epsilon`}{\mjseqn{\lbrace\,\epsilon\,\rbrace}
#'       Standard deviation of the standard normally distributed
#'       phase shift in the seasonally forced transmission rate.
#'     }
#'   }
#'
#' The helper function [make_par_list()] can be used to construct
#' a list conforming to these default requirements.
#'
#' @param par_list A list of parameter values with numeric scalar elements:
#'
#'   \describe{
#'     \item{`dt_days`}{\mjseqn{\lbrace\,\Delta t\,\rbrace}
#'       Observation interval in days.
#'     }
#'     \item{`N0`}{\mjseqn{\lbrace\,N_0\,\rbrace}
#'       Population size at time \mjseqn{t = 0} years.
#'     }
#'     \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'       Number of susceptible individuals at time \mjseqn{t = 0} years.
#'     }
#'     \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'       Number of exposed (infected but not infectious) individuals
#'       at time \mjseqn{t = 0} years. Necessary only if `model = "seir"`.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Number of infectious individuals at time \mjseqn{t = 0} years.
#'     }
#'     \item{`tlat`}{\mjseqn{\lbrace\,t_\text{lat}/\Delta t\,\rbrace}
#'       Mean latent period of the disease of interest
#'       in units \mjseqn{\Delta t}. If `model = "sir"`,
#'       then the removal rate \mjseqn{\gamma} is set equal to
#'       \mjseqn{1 / (t_\text{lat} + t_\text{inf})} (see Details 1).
#'     }
#'     \item{`tinf`}{\mjseqn{\lbrace\,t_\text{inf}/\Delta t\,\rbrace}
#'       Mean infectious period of the disease of interest
#'       in units \mjseqn{\Delta t}. If `model = "sir"`,
#'       then the removal rate \mjseqn{\gamma} is set equal to
#'       \mjseqn{1 / (t_\text{lat} + t_\text{inf})} (see Details 1).
#'     }
#'     \item{`prep`}{\mjseqn{\lbrace\,p_\text{rep}\,\rbrace}
#'       Probability that an infection is reported.
#'     }
#'   }
#'
#'   Additional elements must be included as necessary to ensure that
#'   the functions `nu`, `mu`, and `beta` (see below) can be evaluated.
#' @param n An integer scalar. The number of observations in the
#'   simulated time series, which will have time points
#'   \mjseqn{t_i = i \Delta t} (for \mjseqn{i = 0,\ldots,n-1}).
#' @param with_ds A logical scalar. Should the simulation include
#'   demographic stochasticity? The simulation is performed using
#'   [adaptivetau::ssa.adaptivetau()] if `TRUE` and [deSolve::ode()]
#'   if `FALSE`.
#' @param model A character scalar, either `"sir"` or `"seir"`,
#'   indicating a system of equations to simulate (see Details 1).
#' @param mu A function taking as input a numeric vector `s`
#'   and a list `par_list` and returning as output a numeric vector
#'   of the same length as `s` containing the value of
#'   \mjseqn{f(s) = \mu(s \Delta t) \Delta t} at each element of `s`.
#'   The function must define \mjseqn{f(s)} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\mu(t) \Delta t}
#'   is the per capita natural mortality rate at time \mjseqn{t}
#'   expressed per unit \mjseqn{\Delta t}.
#' @param nu A function taking as input a numeric vector `s`
#'   and a list `par_list` and returning as output a numeric vector
#'   of the same length as `s` containing the value of
#'   \mjseqn{f(s) = \nu(s \Delta t) \Delta t} at each element of `s`.
#'   The function must define \mjseqn{f(s)} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\nu(t) \Delta t}
#'   is the birth rate (relative to \mjseqn{N_0}) at time \mjseqn{t}
#'   expressed per unit \mjseqn{\Delta t}.
#' @param beta A function taking as input a numeric vector `s`
#'   and a list `par_list` and returning as output a numeric vector
#'   of the same length as `s` containing the value of
#'   \mjseqn{f(s) = \beta(s \Delta t) \Delta t} at each element of `s`.
#'   The function must define \mjseqn{f(s)} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\beta(t) \Delta t}
#'   is the transmission rate at time \mjseqn{t}
#'   expressed per unit \mjseqn{\Delta t} per susceptible individual
#'   per infectious individual.
#' @param delay_dist A numeric vector. The distribution of the integer
#'   number of observation intervals between infection and reporting.
#'   `delay_dist[i]` is the probability that an infection that is
#'   eventually reported is reported after `i - 1` observation intervals.
#'   `delay_dist` is replaced with `delay_dist / sum(delay_dist)` in the
#'   event that `sum(delay_dist) != 1`.
#'
#' @return
#' A data frame with `n` rows corresponding to equally spaced time points
#' \mjseqn{t_i = i \Delta t} (for \mjseqn{i = 0, \ldots, n-1}), and numeric
#' columns:
#'
#' \describe{
#'   \item{`t`}{\mjseqn{\lbrace\,t_i / \Delta t\,\rbrace}
#'     Time in units \mjseqn{\Delta t}. Equal to `0:(n - 1)`.
#'   }
#'   \item{`t_years`}{\mjseqn{\lbrace\,t_i\,\rbrace}
#'     Time in years. Equal to `(0:(n - 1)) * par_list$dt_days * / 365`.
#'   }
#'   \item{`mu`}{\mjseqn{\lbrace\,\mu(t_i) \Delta t\,\rbrace}
#'     Per capita natural mortality rate expressed
#'     per unit \mjseqn{\Delta t}. Equal to `mu(0:(n - 1))`.
#'   }
#'   \item{`nu`}{\mjseqn{\lbrace\,\nu(t_i) \Delta t\,\rbrace}
#'     Birth rate (relative to \mjseqn{N_0}) expressed
#'     per unit \mjseqn{\Delta t}. Equal to `nu(0:(n - 1))`.
#'   }
#'   \item{`beta`}{\mjseqn{\lbrace\,\beta(t_i) \Delta t\,\rbrace}
#'     Transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'     Equal to `beta(0:(n - 1))`.
#'   }
#'   \item{`N`}{\mjseqn{\lbrace\,N(t_i)\,\rbrace}
#'     Population size.
#'   }
#'   \item{`S`}{\mjseqn{\lbrace\,S(t_i)\,\rbrace}
#'     Number of susceptible individuals.
#'   }
#'   \item{`E`}{\mjseqn{\lbrace\,E(t_i)\,\rbrace}
#'     Number of exposed (infected but not infectious) individuals.
#'     Included only if `model = "seir"`.
#'   }
#'   \item{`I`}{\mjseqn{\lbrace\,I(t_i)\,\rbrace}
#'     Number of infectious individuals.
#'   }
#'   \item{`Z`}{\mjseqn{\lbrace\,Z(t_i)\,\rbrace}
#'     Incidence. `Z[i]` is the number of infections
#'     between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`Zrep`}{\mjseqn{\lbrace\,Z_\text{rep}(t_i)\,\rbrace}
#'     Incidence conditional on reporting. `Zrep[i]` is the
#'     number of infections between times `t[i-1]` and `t[i]`
#'     that are eventually reported.
#'   }
#'   \item{`C`}{\mjseqn{\lbrace\,C(t_i)\,\rbrace}
#'     Reported incidence. `C[i]` is the number of infections reported
#'     between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`B`}{\mjseqn{\lbrace\,B(t_i)\,\rbrace}
#'     Births. `B[i]` is the number of births
#'     between times `t[i-1]` and `t[i]`.
#'   }
#' }
#'
#' The data frame has attributes `call` and `arg_list`, making it
#' reproducible with `eval(call)` or `do.call(make_data, arg_list)`,
#' possibly preceded by a call to [base::set.seed()].
#'
#' @examples
#' # Stochastic simulation of SEIR model
#' par_list <- make_par_list(
#'   epsilon = pi / 4, # environmental stochasticity
#'   prep = 0.5, # random under-reporting of infections
#'   model = "seir"
#' )
#' set.seed(1422)
#' df <- make_data(par_list,
#'   with_ds = TRUE, # demographic stochasticity
#'   model = "seir",
#'   delay_dist = dnbinom(0:10, mu = 2, size = 4) # random delays in reporting
#' )
#' head(df)
#'
#' # Deterministic simulation of SIR model
#' par_list <- make_par_list(
#'   epsilon = 0,
#'   prep = 1,
#'   model = "sir"
#' )
#' df <- make_data(par_list,
#'   with_ds = FALSE,
#'   model = "sir",
#'   delay_dist = c(1)
#' )
#' head(df)
#'
#' @seealso [make_par_list()], [make_beta()]
#' @export
#' @importFrom stats rbinom
#' @importFrom deSolve ode
#' @importFrom adaptivetau ssa.adaptivetau
make_data <- function(par_list = make_par_list(model = "sir"),
                      n = 1e03,
                      with_ds = FALSE,
                      model = "sir",
                      mu = function(s, par_list) rep(par_list$muconst, length(s)),
                      nu = mu,
                      beta = make_beta(par_list$epsilon, n),
                      delay_dist = c(1)) {
  ## Setup ---------------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Time points
  n <- floor(n)
  times <- 0:(n - 1)

  # Some derived quantities
  if (model == "sir") {
    gamma <- 1 / (par_list$tlat + par_list$tinf)
  } else if (model == "seir") {
    sigma <- 1 / par_list$tlat
    gamma <- 1 / par_list$tinf
  } else {
    stop("`model` must be one of `\"sir\"` or `\"seir\"`.")
  }

  # Normalize `delay_dist`
  delay_dist <- delay_dist / sum(delay_dist)


  ## Simulate S(E)IR equations ... ---------------------------------------
  ## if with demographic stochasticity

  if (with_ds) {

    if (model == "sir") {

      # Initial state
      x_init <- with(par_list, {
        c(
          S = floor(S0),
          I = floor(I0),
          R = floor(N0) - floor(S0) - floor(I0),
          Bcum = 0,
          Zcum = 0
        )
      })

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
        m <- mu(t, params)
        c(
          nu(t, params) * params$N0,         # birth
          beta(t, params) * x["S"] * x["I"], # infection
          gamma * x["I"] * (x["I"] > 1),     # removal
          m * x["S"],                        # natural mortality
          m * x["I"] * (x["I"] > 1),
          m * x["R"]
        )
      }

    } else if (model == "seir") {

      # Initial state
      x_init <- with(par_list, {
        c(
          S = floor(S0),
          E = floor(E0),
          I = floor(I0),
          R = floor(N0) - floor(S0) - floor(E0) - floor(I0),
          Bcum = 0,
          Zcum = 0
        )
      })

      # Transition events
      event_list <- list(
        c(S = 1, Bcum = 1),         # birth
        c(S = -1, E = 1, Zcum = 1), # infection
        c(E = -1, I = 1),           # onset of infectiousness
        c(I = -1, R = 1),           # removal
        c(S = -1),                  # natural mortality
        c(E = -1),
        c(I = -1),
        c(R = -1)
      )

      # Transition event rates
      compute_event_rates <- function(x, params, t) {
        m <- mu(t, params)
        c(
          nu(t, params) * params$N0,              # birth
          beta(t, params) * x["S"] * x["I"],      # infection
          sigma * x["E"],                         # onset of infectiousness
          gamma * x["I"] * (x["E"] + x["I"] > 1), # removal
          m * x["S"],                             # natural mortality
          m * x["E"] * (x["E"] + x["I"] > 1),
          m * x["I"] * (x["E"] + x["I"] > 1),
          m * x["R"]
        )
      }

    }

    # Generate a realization of the stochastic process
    df <- as.data.frame(
      ssa.adaptivetau(
        init.values = x_init,
        transitions = event_list,
        rateFunc    = compute_event_rates,
        params      = par_list,
        tf          = n - 1, # final time point
        tl.params   = list(
          epsilon     = 0.05,
          delta       = 0.05,
          maxtau      = 0.5, # adaptive time step must not exceed 1
          extraChecks = TRUE
        )
      )
    )
    colnames(df)[1] <- "t"

    ## NOTE: `ssa.adaptivetau()` returns time series with unequal spacing
    ##       (time step varies between 0 and `maxtau`), but we desire
    ##       equal spacing (time step equal to 1 as in `times`)

    # For each time in `times`, find the index of the last previous time
    # in `df$t`
    ind_state_out <- sapply(times, function(s) max(which(df$t <= s)))

    # States at those times are states at times `times`
    df <- df[ind_state_out, ]
    df$t <- times
    rownames(df) <- NULL


  ## Simulate S(E)IR equations ... ---------------------------------------
  ## if without demographic stochasticity

  } else {
    if (model == "sir") {

      # Initial state
      x_init <- with(par_list, {
        c(
          S    = S0,
          logI = log(I0),
          R    = N0 - S0 - I0,
          Bcum = 0,
          Zcum = 0
        )
      })

      # System of SIR equations
      compute_ode_rates <- function(t, y, parms) {
        m <- mu(t, parms)
        b <- beta(t, parms)
        y["I"] <- exp(y["logI"])
        dBcum <- nu(t, parms) * parms$N0
        dZcum <- b * y["S"] * y["I"]
        dS <- dBcum - dZcum - m * y["S"]
        dlogI <- b * y["S"] - gamma - m
        dR <- gamma * y["I"] - m * y["R"]
        list(c(dS, dlogI, dR, dBcum, dZcum))
      }

    } else if (model == "seir") {

      # Initial state
      x_init <- with(par_list, {
        c(
          S    = S0,
          logE = log(E0),
          logI = log(I0),
          R    = N0 - S0 - E0 - I0,
          Bcum = 0,
          Zcum = 0
        )
      })

      # System of SEIR equations
      compute_ode_rates <- function(t, y, parms) {
        m <- mu(t, parms)
        b <- beta(t, parms)
        y["E"] <- exp(y["logE"])
        y["I"] <- exp(y["logI"])
        dBcum <- nu(t, parms) * parms$N0
        dZcum <- b * y["S"] * y["I"]
        dS <- dBcum - dZcum - m * y["S"]
        dlogE <- b * y["S"] * y["I"] / y["E"] - sigma - m
        dlogI <- sigma * y["E"] / y["I"] - gamma - m
        dR <- gamma * y["I"] - m * y["R"]
        list(c(dS, dlogE, dlogI, dR, dBcum, dZcum))
      }

    }

    # Numerically integrate the system of S(E)IR equations
    df <- as.data.frame(
      ode(
        y     = x_init,
        times = times,
        func  = compute_ode_rates,
        parms = par_list
      )
    )
    colnames(df)[1] <- "t"
    if (model == "sir") {
      df <- transform(df, I = exp(logI))
    } else if (model == "seir") {
      df <- transform(df, E = exp(logE), I = exp(logI))
    }

  }


  ## Append other desired variables --------------------------------------

  df <- transform(df,
    t_years  = t * par_list$dt_days / 365,
    mu       = mu(t, par_list),
    nu       = nu(t, par_list),
    beta     = beta(t, par_list),
    N        = S + (if (model == "seir") E else 0) + I + R,
    Z        = c(NA, diff(Zcum)), # `Z` from `Zcum` by first differences
    B        = c(NA, diff(Bcum))  # `B` from `Bcum` by first differences
  )


  ## Introduce observation error -----------------------------------------

  # `Zrep` from `Z` by binomial sampling
  df$Zrep <- c(NA,
    stats::rbinom(
      n    = nrow(df) - 1,    # number of experiments
      size = round(df$Z[-1]), # numbers of Bernoulli trials
      p    = par_list$prep    # success probability
    )
  )

  # `C` from `Zrep` by binning forward in time
  convol_out <- convol(df$Zrep[-1], delay_dist, n = 1)
  df$C <- c(NA, convol_out$simulation[1:(nrow(df)-1), 1])

  if (model == "sir") {
    df <- df[, c("t", "t_years", "mu", "nu", "beta",
                 "N", "S", "I", "Z", "Zrep", "C", "B")]
  } else if (model == "seir") {
    df <- df[, c("t", "t_years", "mu", "nu", "beta",
                 "N", "S", "E", "I", "Z", "Zrep", "C", "B")]
  }
  structure(df, call = match.call(), arg_list = arg_list)
}
