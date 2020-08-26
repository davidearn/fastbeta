#' \loadmathjax
#' Simulate epidemic time series data
#'
#' @description
#' Simulates epidemic time series data using a system of S(E)IR equations
#' and a supplied parametrization (see Details 1). Observations are recorded
#' at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}. Users can
#' specify any time-varying birth, death, and transmission rates and any
#' discrete distribution of the time from infection to reporting. By default,
#' the vital rates are constant and the transmission rate is seasonally forced.
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
#' \mjseqn{\gamma = 1 / t_\text{inf}}.
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
#'       Number of exposed individuals at time \mjseqn{t = 0} years.
#'       Necessary only if `model = "seir"`.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Number of infected (`model = "sir"`) or infectious (`model = "seir"`)
#'       individuals at time \mjseqn{t = 0} years.
#'     }
#'     \item{`tlat`}{\mjseqn{\lbrace\,t_\text{lat}/\Delta t\,\rbrace}
#'       Mean latent period of the disease of interest
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`tinf`}{\mjseqn{\lbrace\,t_\text{inf}/\Delta t\,\rbrace}
#'       Mean infectious period of the disease of interest
#'       in units \mjseqn{\Delta t}.
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
#'   expressed per unit \mjseqn{\Delta t}.
#' @param delay_dist A numeric vector. The distribution of the
#'   number of observations intervals between infection and reporting.
#'   `delay_dist[i]` is the probability that an infection that is
#'   eventually reported is reported after `i - 1` observation intervals.
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
#'     Number of exposed individuals. Included only if `model = "seir"`.
#'   }
#'   \item{`I`}{\mjseqn{\lbrace\,I(t_i)\,\rbrace}
#'     Number of infected (`model = "sir"`) or infectious (`model = "seir"`)
#'     individuals.
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
#' preceded by a call to `set.seed()`.
#'
#' @examples
#' # Stochastic simulation of SIR model
#' par_list <- make_par_list(
#'   epsilon = pi / 4, # environmental stochasticity
#'   prep = 0.5, # random under-reporting of infections
#'   model = "seir"
#' )
#' set.seed(1422)
#' df <- make_data(par_list,
#'   with_ds = TRUE, # demographic stochasticity
#'   model = "sir",
#'   delay_dist = dnbinom(0:10, mu = 2, size = 4) /
#'     pnbinom(10, mu = 2, size = 4) # random delays in reporting
#' )
#' head(df)
#'
#' # Deterministic simulation of SEIR model
#' par_list <- make_par_list(
#'   epsilon = 0,
#'   prep = 1,
#'   model = "seir"
#' )
#' df <- make_data(par_list,
#'   with_ds = FALSE,
#'   model = "seir",
#'   delay_dist = c(1)
#' )
#' head(df)
#'
#' @seealso [make_par_list()], [make_beta()]
#' @export
#' @importFrom stats rbinom
#' @importFrom deSolve ode
#' @importFrom adaptivetau ssa.adaptivetau
make_data <- function(par_list = make_par_list(),
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
    Z        = c(NA, diff(Zcum)), # `Z` from `Zcum` via first differences
    B        = c(NA, diff(Bcum))  # `B` from `Bcum` via first differences
  )


  ## Introduce observation error -----------------------------------------

  # `Zrep` from `Z` via binomial sampling
  df$Zrep <- c(NA,
    stats::rbinom(
      n    = nrow(df) - 1,    # number of experiments
      size = round(df$Z[-1]), # numbers of Bernoulli trials
      p    = par_list$prep    # success probability
    )
  )

  # `C` from `Zrep` via binning forward in time
  counts <- df$Zrep[-1]
  delays <- sample(
    x       = 0:(length(delay_dist) - 1),
    size    = sum(counts),
    replace = TRUE,
    prob    = delay_dist
  )
  index_counts <- rep(seq_along(counts), times = counts)
  index_counts_with_delays <- index_counts + delays
  counts_with_delays <- rep(0, max(index_counts_with_delays))
  counts_with_delays[sort(unique(index_counts_with_delays))] <-
    table(index_counts_with_delays)
  df$C <- c(NA, counts_with_delays[1:(nrow(df)-1)])

  if (model == "sir") {
    df <- df[, c("t", "t_years", "mu", "nu", "beta",
                 "N", "S", "I", "Z", "Zrep", "C", "B")]
  } else if (model == "seir") {
    df <- df[, c("t", "t_years", "mu", "nu", "beta",
                 "N", "S", "E", "I", "Z", "Zrep", "C", "B")]
  }
  structure(df, call = match.call(), arg_list = arg_list)
}

#' \loadmathjax
#' Define a seasonally forcing function for simulations
#'
#' @description
#' Defines a seasonal forcing function with environmental noise
#' that can be assigned to the argument `beta` of [make_data()].
#' 
#' @details
#' A noise process \mjseqn{\lbrace(i;X_i)\rbrace_{i=0}^{n-1}}
#' with \mjseqn{X_i \sim \mathrm{Normal}(0,\epsilon^2)} is
#' realized and linearly interpolated, yielding a randomly generated
#' function \mjseqn{a(s)} defined continuously on the interval
#' \mjseqn{\lbrack 0,n-1 \rbrack}. A seasonal forcing function
#' `f()` is then defined according to
#'
#' \mjsdeqn{f(s) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi s \Delta t}{\text{365 days}} + a(s)\unicode{x1D7D9}_{\lbrack 0,n-1 \rbrack}\right) \right\rbrack \Delta t}
#'
#' and returned as output. `f()` takes as arguments:
#'
#' \describe{
#'   \item{`s`}{A numeric vector listing values of \mjseqn{s}
#'     at which to evaluate \mjseqn{f(s)}.
#'   }
#'   \item{`par_list`}{A list with numeric scalar elements
#'     `dt_days`, `beta_mean`, and `alpha` listing  values for
#'     \mjseqn{\Delta t}, \mjseqn{\langle\beta\rangle}, and
#'     \mjseqn{\alpha}, respectively.
#'   }
#' }
#'
#' `par_list` does not need an element `epsilon` (giving a value
#' for \mjseqn{\epsilon}), because random number generation and
#' linear interpolation are performed in the enclosing environment
#' of `f()` (the execution environment of `make_beta()`) using the
#' value of `epsilon` found there. That is, `f()` does not need to
#' know `epsilon`, because it does no sampling or interpolating of
#' its own and instead uses the definition of the interpolant
#' \mjseqn{a(s)} that it finds in its enclosing environment.
#' The upshot is that the function `f()` is itself randomly generated
#' (reproducible with [base::set.seed()]),
#' whereas the output of `f()` is determined
#' (reproducible without [base::set.seed()]).
#'
#' @param epsilon A numeric scalar. The standard deviation of
#'   the noise process (see Details).
#' 
#' @param n An integer scalar. The number of observations in
#'   the noise process (see Details).
#' 
#' @return
#' A function with arguments `s` and `par_list` (see Details).
#'
#' @examples
#' set.seed(1734)
#' epsilon <- pi / 4
#' n <- 1e03
#' f <- make_beta(epsilon, n)
#' formals(f)
#' body(f)
#' names(as.list(environment(f)))
#' a <- get("a", envir = environment(f))
#' s <- 0:(n - 1)
#' plot(s, a(s), type = "l", col = "blue")
#'
#' @seealso [make_data()]
#' @export
#' @importFrom stats approxfun rnorm
make_beta <- function(epsilon, n) {
  n <- floor(n)
  a <- approxfun(
    x      = 0:(n - 1),
    y      = rnorm(n, mean = 0, sd = epsilon),
    method = "linear"
  )
  function(s, par_list) {
    a_val <- ifelse(s < 0 | s > n - 1, 0, a(s))
    one_year <- 365 / par_list$dt_days
    par_list$beta_mean *
      (1 + par_list$alpha * cos(2 * pi * s / one_year + a_val))
  }
}

#' \loadmathjax
#' Create a list of parameter values for simulations
#'
#' @description
#' Creates a list of parameter values that can be assigned
#' to the argument `par_list` of [make_data()] if defaults
#' are accepted for its arguments `mu`, `nu`, and `beta`.
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
#' or the identity 
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
#' Choosing \mjseqn{n} such that \mjseqn{n \Delta t \sim 1000} years is
#' typically enough to ensure that \mjseqn{\big(S(0),I(0),E(0)\big)} is
#' near the attractor of the system, which can be desirable when simulating
#' epidemic time series data (e.g., using [make_data()]).
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
#'   A numeric scalar. Birth rate (relative to \mjseqn{N_0})
#'   expressed per unit \mjseqn{\Delta t} *and* per capita
#'   natural mortality rate expressed per unit \mjseqn{\Delta t}.
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
#' @param prep \mjseqn{\lbrace\,p_\text{rep}\,\rbrace}
#'   A numeric scalar. Probability that an infection is reported.
#' @param model A character scalar, either `"sir"` or `"seir"`,
#'   indicating a system of equations to be numerically integrated
#'   (see Details 2).
#' @param n \mjseqn{\lbrace\,n\,\rbrace}
#'   An integer scalar (positive). A system of SIR (`model = "sir"`)
#'   or SEIR (`model = "seir"`) equations will be numerically integrated
#'   between times \mjseqn{t = -(n-1) \Delta t} and \mjseqn{t = 0} years
#'   to obtain values for \mjseqn{S_0}, \mjseqn{E_0} (`model = "seir"`),
#'   and \mjseqn{I_0} (see Details 2 and 3).
#'
#' @return
#' A list of the arguments in the function call,
#' excluding `model` and `n` but including these additional elements:
#'
#' \describe{
#'   \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'     Mean (long-term average) of the seasonally forced
#'     transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'   }
#'   \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'     Number of susceptible individuals at time \mjseqn{t = 0} years.
#'   }
#'   \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'     Number of exposed individuals at time \mjseqn{t = 0} years.
#'     Included only if `model = "seir"`.
#'   }
#'   \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'     Number of infected (`model = "sir"`) or infectious (`model = "seir"`)
#'     individuals at time \mjseqn{t = 0} years.
#'   }
#' }
#'
#' This list can be assigned to the argument `par_list` of `[make_data()]`
#' if defaults are accepted for its arguments `mu`, `nu`, and `beta`.
#' In this case, care must be taken to ensure that the argument `model`
#' of [make_data()] matches that of `make_par_list()`.
#' 
#' @examples
#' # Creates a list for measles in England by default
#' par_list <- make_par_list()
#' unlist(par_list)
#'
#' # Time (rate) parameters should be specified
#' # in units (per unit) of the observation interval
#' dt_days <- 7
#' tlat_days <- 5
#' tinf_weeks <- 8 / 7
#' muconst_peryear <- 0.04
#' par_list <- make_par_list(
#'   dt_days = dt_days,
#'   tlat    = tlat_days / dt_days,
#'   tinf    = tinf_weeks * 7 / dt_days,
#'   muconst = muconst_peryear * dt_days / 365
#' )
#' unlist(par_list)
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
                          prep     = 1,
                          model    = "seir",
                          n        = 1000 * 365 / dt_days) {
  # Some derived quantities
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

  # Time points
  n <- floor(n)
  times <- -(n - 1):0

  if (model == "sir") {
    
    # Initial state
    x_init <- c(
      S    = N0 * (1 / Rnaught),
      logI = log(N0 * (1 - 1 / Rnaught) * (muconst / (gamma + muconst)))
    )

    # System of SIR equations
    compute_ode_rates <- function(t, y, parms) {
      y["I"] <- exp(y["logI"])
      beta <- beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
      dS <- muconst * N0 - beta * y["S"] * y["I"] - muconst * y["S"]
      dlogI <- beta * y["S"] - gamma - muconst
      list(c(dS, dlogI))
    }
    
  } else if (model == "seir") {

    # Initial state
    x_init <- c(
      S    = N0*(1/Rnaught),
      logE = log(N0*(1-1/Rnaught)*(muconst/(sigma+muconst))),
      logI = log(N0*(1-1/Rnaught)*(sigma/(sigma+muconst))*(muconst/(gamma+muconst)))
    )

    # System of SEIR equations
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

  # Numerically integrate the system of S(E)IR equations
  df <- as.data.frame(
    ode(
      y     = x_init,
      times = times,
      func  = compute_ode_rates,
      parms = NULL, # found in enclosing environment of `compute_seir_rates()`
      hmax  = 1 # avoids error when `length(times) = 1`
    )
  )

  # Assign final values of `S`, `E`, and `I`
  S0 <- df[n, "S"]
  if (model == "seir") {
    E0 <- exp(df[n, "logE"])
  }
  I0 <- exp(df[n, "logI"])

  if (model == "sir") {
    as.list(environment())[c("dt_days", "N0", "S0", "I0",
                             "tlat", "tinf", "muconst", "Rnaught",
                             "beta_mean", "alpha", "epsilon", "prep")]
  } else if (model == "seir") {
    as.list(environment())[c("dt_days", "N0", "S0", "E0", "I0",
                             "tlat", "tinf", "muconst", "Rnaught",
                             "beta_mean", "alpha", "epsilon", "prep")]
  }
}
