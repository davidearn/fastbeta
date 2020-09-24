#' \loadmathjax
#' Simulate epidemic time series data
#'
#' @description
#' Simulates epidemic time series data using a system of S(E)IR equations
#' (see Details 1). Observations are recorded at equally spaced time points
#' \mjseqn{t_i = i \Delta t}. Users can specify any time-varying rates of
#' birth, death, and transmission, any time-varying probability that an
#' infection is reported, and any discrete distribution of the time from
#' infection to reporting. By default, the vital rates are constant, the
#' transmission rate is seasonally forced, and infections are reported
#' with a constant probability and zero delay.
#'
#' @details
#' # Details
#'
#' ## 1. Simulation model
#'
#' `make_data()` simulates epidemic time series data using a system
#' of SIR or SEIR equations that includes additional equations for
#' cumulative births and cumulative incidence. If `model = "sir"`,
#' then the system is
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu(t)\,, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta(t) S I\,, \end{align*}}
#'
#' where \mjseqn{\gamma = 1 / (t_\text{lat} + t_\text{inf})}.
#' If `model = "seir"`, then the system is
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}E}{\text{d}t} &= \beta(t) S I - \sigma E - \mu(t) E\,, \cr \frac{\text{d}I}{\text{d}t} &= \sigma E - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu(t)\,, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta(t) S I\,, \end{align*}}
#'
#' where \mjseqn{\sigma = 1 / t_\text{lat}} and
#' \mjseqn{\gamma = 1 / t_\text{inf}}.
#'
#' `make_data()` generates observations of the system at equally spaced
#' time points \mjseqn{t_i = i \Delta t} (\mjseqn{i = 0,\ldots,n-1}) by
#' either
#' (i) numerically integrating the system of ODE using [deSolve::ode()]
#' (`with_ds = FALSE`) or
#' (ii) realizing a corresponding continuous-time stochastic process
#' using [adaptivetau::ssa.adaptivetau()] (`with_ds = TRUE`).
#' Both methods use initial state
#'
#' \mjsdeqn{\begin{bmatrix} S(0) \cr I(0) \cr R(0) \cr B_\text{cum}(0) \cr Z_\text{cum}(0) \end{bmatrix} = \begin{bmatrix} S_0 \cr I_0 \cr N_0 - S_0 - I_0 \cr 0 \cr 0 \end{bmatrix}}
#'
#' for `model = "sir"` and
#'
#' \mjsdeqn{\begin{bmatrix} S(0) \cr E(0) \cr I(0) \cr R(0) \cr B_\text{cum}(0) \cr Z_\text{cum}(0) \end{bmatrix} = \begin{bmatrix} S_0 \cr E_0 \cr I_0 \cr N_0 - S_0 - E_0 - I_0 \cr 0 \cr 0 \end{bmatrix}}
#'
#' for `model = "seir"`.
#'
#' The latter method (`with_ds = TRUE`) defines event probabilities as
#' proportional to terms in the ODE, modeling **demographic stochasticity**.
#' Disease fadeout in simulations with demographic stochasticity is prevented
#' by setting the probability of a decrease in the infected population size
#' to zero whenever the infected population size is one.
#'
#' `make_data()` calculates births \mjseqn{B} and incidence \mjseqn{Z}
#' from cumulative births \mjseqn{B_\text{cum}} and cumulative incidence
#' \mjseqn{Z_\text{cum}} by first differences:
#'
#' \mjsdeqn{\begin{align*} B(t_i) &= B_\text{cum}(t_i) - B_\text{cum}(t_{i-1}), \cr Z(t_i) &= Z_\text{cum}(t_i) - Z_\text{cum}(t_{i-1}). \end{align*}}
#'
#' This calculation leaves \mjseqn{B(t_0)} and \mjseqn{Z(t_0)} undefined.
#' For simplicity, they are assigned the values of \mjseqn{B(t_1)} and
#' \mjseqn{Z(t_1)}, respectively.
#'
#' `make_data()` then does binomial sampling to generate, for each
#' observation of \mjseqn{Z}, a number of infections that is eventually
#' reported:
#'
#' \mjsdeqn{Z_\text{rep}(t_i) \sim \mathrm{Binomial}\big(\mathrm{nint}(Z(t_i)),p_i\big)\,.}
#'
#' Here, \mjseqn{\mathrm{nint}(Z(t_i))} is \mjseqn{Z(t_i)} rounded to the
#' nearest integer and \mjseqn{p_i} is the probability that an infection
#' between times \mjseqn{t_{i-1}} and \mjseqn{t_i} (counted by \mjseqn{Z(t_i)})
#' is eventually reported. Finally, it generates reported incidence \mjseqn{C}
#' by sampling from a supplied reporting delay distribution and binning
#' infections counted by \mjseqn{Z_\text{rep}} forward in time. If the
#' maximum possible delay is \mjseqn{b \Delta t} and \mjseqn{b > 0},
#' then, for simplicity, `make_data()` simulates
#'
#' \mjsdeqn{Z_\text{rep}(t_i) \sim \mathrm{Binomial}\big(\mathrm{nint}(Z(t_0)),p_0\big)\,.}
#'
#' for \mjseqn{i = -b,\ldots,-1}, in order to ensure that \mjseqn{C} is based
#' on complete information about \mjseqn{Z_\text{rep}}. Compared to \mjseqn{Z},
#' \mjseqn{C} models **observation error**, i.e., under-reporting of infections
#'  with delays between infection and reporting.
#'
#' ## 2. Default parametrization
#'
#' By default, \mjseqn{\mu} and \mjseqn{\nu} are modeled as constant
#' and equal to \mjseqn{\mu_\text{c}} and \mjseqn{\mu_\text{c} N_0},
#' respectively, while \mjseqn{\beta} is modeled by a seasonal forcing
#' function with environmental noise. Specifically,
#'
#' \mjsdeqn{\beta(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{365 days}} + \phi(t;\epsilon^2)\right) \right\rbrack\,,}
#'
#' where \mjseqn{\phi(t;\epsilon^2)} is the linear interpolant of noise
#' \mjseqn{\lbrace(t_i;\Phi_i)\rbrace_{i=0}^{n-1}} with
#'
#' \mjsdeqn{\Phi_i \sim \mathrm{Normal}(0,\epsilon^2)\,.}
#'
#' Care is taken to ensure that random number generation occurs in the
#' enclosing environment rather than the body of the function assigned
#' by default to argument `beta`. This ensures that the function has
#' noisy but non-random output, which is helpful for reproducibility.
#' For details, see helper function [make_beta()].
#'
#' By default, the reporting probability \mjseqn{p} is constant and equal
#' to \mjseqn{p_\text{c}}, while the reporting delay is fixed equal to zero,
#' so that all infections during a given observation interval that are
#' eventually reported, are reported during the same observation interval,
#' i.e., \mjseqn{C(t_i) = Z_\text{rep}(t_i)} for all \mjseqn{i}.
#'
#' The default specification of `mu`, `nu`, `beta`, and `p` requires
#' that `par_list` contains these additional numeric scalar elements:
#'
#' \describe{
#'   \item{`muconst`}{\mjseqn{\lbrace\,\mu_\text{c} \Delta t\,\rbrace}
#'     Per capita natural mortality rate
#'     expressed per unit \mjseqn{\Delta t}.
#'   }
#'   \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'     Mean (long-term average) of the seasonally forced
#'     transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'   }
#'   \item{`alpha`}{\mjseqn{\lbrace\,\alpha\,\rbrace}
#'     Amplitude of the seasonally forced transmission rate
#'     relative to the mean.
#'   }
#'   \item{`epsilon`}{\mjseqn{\lbrace\,\epsilon\,\rbrace}
#'     Standard deviation of the standard normally distributed
#'     phase shift in the seasonally forced transmission rate.
#'   }
#'   \item{`pconst`}{\mjseqn{\lbrace\,p_\text{c}\,\rbrace}
#'     Probability that an infection is eventually reported.
#'   }
#' }
#'
#' Helper function [make_par_list()] can be used to construct
#' a list conforming to this requirement.
#'
#' @param par_list A list of parameter values with numeric scalar elements:
#'
#'   \describe{
#'     \item{`dt_days`}{\mjseqn{\lbrace\,\Delta t\,\rbrace}
#'       Observation interval in days.
#'     }
#'     \item{`tlat`}{\mjseqn{\lbrace\,t_\text{lat}/\Delta t\,\rbrace}
#'       Mean latent period of the disease of interest
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`tinf`}{\mjseqn{\lbrace\,t_\text{inf}/\Delta t\,\rbrace}
#'       Mean infectious period of the disease of interest
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`N0`}{\mjseqn{\lbrace\,N_0\,\rbrace}
#'       Population size at time \mjseqn{t_0 = 0 \Delta t}.
#'     }
#'     \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'       Number of susceptible individuals at time \mjseqn{t_0 = 0 \Delta t}.
#'     }
#'     \item{`E0`}{\mjseqn{\lbrace\,E_0\,\rbrace}
#'       Number of exposed (infected but not infectious) individuals
#'       at time \mjseqn{t_0 = 0 \Delta t}. Necessary only if `model = "seir"`.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Number of infectious individuals at time \mjseqn{t_0 = 0 \Delta t}.
#'     }
#'   }
#'
#'   Additional elements must be included as necessary to ensure
#'   that arguments `mu`, `nu`, and `beta` are well-defined
#'   (see Details 2).
#' @param n \mjseqn{\lbrace\,n\,\rbrace}
#'   An integer scalar. The number of observations in the
#'   simulated time series, which will have time points
#'   \mjseqn{t_i = i \Delta t} (\mjseqn{i = 0,\ldots,n-1}).
#' @param with_ds A logical scalar. Should the simulation include
#'   demographic stochasticity? The simulation is performed using
#'   [deSolve::ode()] if `FALSE` and [adaptivetau::ssa.adaptivetau()]
#'   if `TRUE`.
#' @param model A character scalar, either `"sir"` or `"seir"`,
#'   indicating a model of disease dynamics (see Details 1).
#' @param mu \mjseqn{\lbrace\,\mu(s \Delta t) \Delta t\,\rbrace}
#'   A function taking as arguments a numeric vector `s` and a list `par_list`
#'   and returning as output a numeric vector of length `length(s)` containing
#'   the value of \mjseqn{\mu(s \Delta t) \Delta t} at each element of `s`.
#'   Must define \mjseqn{\mu(s \Delta t) \Delta t} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\mu(t) \Delta t}
#'   is the per capita natural mortality rate at time \mjseqn{t} expressed
#'   per unit \mjseqn{\Delta t}. Alternatively, a numeric scalar indicating
#'   a constant rate \mjseqn{\mu_\text{c} \Delta t}.
#' @param nu \mjseqn{\lbrace\,\nu(s \Delta t) \Delta t\,\rbrace}
#'   A function taking as arguments a numeric vector `s` and a list `par_list`
#'   and returning as output a numeric vector of length `length(s)` containing
#'   the value of \mjseqn{\nu(s \Delta t) \Delta t} at each element of `s`.
#'   Must define \mjseqn{\nu(s \Delta t) \Delta t} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\nu(t) \Delta t}
#'   is the birth rate at time \mjseqn{t} expressed per unit \mjseqn{\Delta t}.
#'   Alternatively, a numeric scalar indicating a constant rate
#'   \mjseqn{\nu_\text{c} \Delta t}.
#' @param beta \mjseqn{\lbrace\,\beta(s \Delta t) \Delta t\,\rbrace}
#'   A function taking as arguments a numeric vector `s` and a list `par_list`
#'   and returning as output a numeric vector of length `length(s)` containing
#'   the value of \mjseqn{\beta(s \Delta t) \Delta t} at each element of `s`.
#'   Must define \mjseqn{\beta(s \Delta t) \Delta t} for all
#'   \mjseqn{s \in \lbrack 0,n-1 \rbrack}. Here, \mjseqn{\beta(t) \Delta t}
#'   is the transmission rate at time \mjseqn{t} expressed per unit
#'   \mjseqn{\Delta t} per susceptible individual per infectious individual.
#'   Alternatively, a numeric scalar indicating a constant rate
#'   \mjseqn{\beta_\text{c} \Delta t}.
#' @param p \mjseqn{\lbrace\,p_i\,\rbrace}
#'   A numeric vector of length `n` listing the probability \mjseqn{p_i}
#'   (\mjseqn{i = 0,\ldots,n-1}) that an infection between times
#'   \mjseqn{t_{i-1}} and \mjseqn{t_i} is eventually reported.
#'   Alternatively, a numeric scalar indicating a constant probability
#'   \mjseqn{p_\text{c}}.
#' @param delay_dist \mjseqn{\lbrace\,q_i\,\rbrace}
#'   A numeric vector listing the probability \mjseqn{q_i}
#'   (\mjseqn{i = 0,\ldots,b}) that an infection that is eventually
#'   reported is reported after \mjseqn{i} observation intervals.
#'   `delay_dist` is replaced with `delay_dist / sum(delay_dist)`
#'   in the event that `sum(delay_dist) != 1`.
#'
#' @return
#' A data frame with `n` rows corresponding to equally spaced time points
#' \mjseqn{t_i = i \Delta t} (\mjseqn{i = 0,\ldots,n-1}) and 14 numeric
#' variables:
#'
#' \describe{
#'   \item{`t`}{\mjseqn{\lbrace\,t_i / \Delta t\,\rbrace}
#'     Time in units \mjseqn{\Delta t}. Equal to `0:(n-1)`.
#'   }
#'   \item{`t_years`}{\mjseqn{\lbrace\,t_i\,\rbrace}
#'     Time in years. Equal to `(0:(n-1)) * par_list$dt_days * / 365`.
#'   }
#'   \item{`beta`}{\mjseqn{\lbrace\,\beta(t_i) \Delta t\,\rbrace}
#'     Transmission rate expressed per unit \mjseqn{\Delta t}
#'     per susceptible individual per infectious individual.
#'     Equal to `beta(0:(n-1))` if argument `beta` is a function
#'     and `rep(beta, n)` if it is a numeric scalar.
#'   }
#'   \item{`mu`}{\mjseqn{\lbrace\,\mu(t_i) \Delta t\,\rbrace}
#'     Per capita natural mortality rate expressed
#'     per unit \mjseqn{\Delta t}. Equal to `mu(0:(n-1))`
#'     if argument `mu` is a function and `rep(mu, n)`
#'     if it is a numeric scalar.
#'   }
#'   \item{`nu`}{\mjseqn{\lbrace\,\nu(t_i) \Delta t\,\rbrace}
#'     Birth rate expressed per unit \mjseqn{\Delta t}.
#'     Equal to `nu(0:(n-1))` if argument `nu` is a function
#'     and `rep(nu, n)` if it is a numeric scalar.
#'   }
#'   \item{`B`}{\mjseqn{\lbrace\,B(t_i)\,\rbrace}
#'     Births. `B[i]` is the number of births
#'     between times `t[i-1]` and `t[i]`.
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
#'   \item{`p`}{\mjseqn{\lbrace\,p_i\,\rbrace}
#'     Reporting probability. `p[i]` is the probability
#'     that an infection between times `t[i-1]` and `t[i]`
#'     is eventually reported. Equal to the value of `p`
#'     in the function call.
#'   }
#'   \item{`Zrep`}{\mjseqn{\lbrace\,Z_\text{rep}(t_i)\,\rbrace}
#'     Incidence conditional on reporting. `Zrep[i]` is the
#'     number of infections between times `t[i-1]` and `t[i]`
#'     that are eventually reported.
#'   }
#'   \item{`C`}{\mjseqn{\lbrace\,C(t_i)\,\rbrace}
#'     Reported incidence. `C[i]` is the number of infections reported
#'     between times `t[i-1]` and `t[i]`.
#'     If `b = max(which(delay_dist > 0)) - 1` and `b > 0`,
#'     then attempts are made to ensure that `C[1:b]` have
#'     the correct scale (see Details 1).
#'   }
#' }
#'
#' The data frame has attributes `call` and `arg_list` and
#' is reproducible with `eval(call)` preceded by a call to
#' [base::set.seed()].
#'
#' @examples
#' ## Stochastic simulation of SEIR model
#' pl <- make_par_list(
#'   epsilon = 0.8, # environmental stochasticity
#'   pconst = 0.4, # random under-reporting of infections
#'   model = "seir"
#' )
#' set.seed(17090019)
#' df <- make_data(pl,
#'   with_ds = TRUE, # demographic stochasticity
#'   model = "seir",
#'   delay_dist = dnbinom(0:10, mu = 2, size = 4) # random delays in reporting
#' )
#' head(df)
#'
#' ## Deterministic simulation of SEIR model
#' pl <- make_par_list(
#'   epsilon = 0,
#'   pconst = 1,
#'   model = "sir"
#' )
#' df <- make_data(pl,
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
make_data <- function(par_list   = make_par_list(model = "sir"),
                      n          = 1000L,
                      with_ds    = FALSE,
                      model      = "sir",
                      mu         = par_list$muconst,
                      nu         = par_list$muconst * par_list$N0,
                      beta       = make_beta(par_list$epsilon, n),
                      p          = par_list$pconst,
                      delay_dist = c(1)) {
  ### Setup ------------------------------------------------------------

  ## Save arguments in a list
  arg_list <- as.list(environment())

  ## Time points
  n <- floor(n)
  times <- 0:(n-1)

  ## Some derived quantities
  if (model == "sir") {
    gamma <- 1 / (par_list$tlat + par_list$tinf)
  } else if (model == "seir") {
    sigma <- 1 / par_list$tlat
    gamma <- 1 / par_list$tinf
  } else {
    stop("`model` must be one of `\"sir\"` or `\"seir\"`.")
  }

  ## Map numeric `mu`, `nu`, and `beta` to constant functions
  if (is.numeric(mu)) {
    muconst <- mu[1]
    mu <- function(s, par_list) rep(muconst, length(s))
  }
  if (is.numeric(nu)) {
    nuconst <- nu[1]
    nu <- function(s, par_list) rep(nuconst, length(s))
  }
  if (is.numeric(beta)) {
    betaconst <- beta[1]
    beta <- function(s, par_list) rep(betaconst, length(s))
  }

  ## Map `p` with length not equal to `n`
  ## to a constant vector of length `n`
  if (length(p) != n) {
    pconst <- p[1]
    p <- rep(pconst, n)
  }

  ## Normalize `delay_dist`
  delay_dist <- delay_dist / sum(delay_dist)


  ### Simulate S(E)IR equations ... ------------------------------------
  ### if with demographic stochasticity

  if (with_ds) {

    if (model == "sir") {

      ## Initial state
      x_init <- with(par_list, {
        c(
          S = floor(S0),
          I = floor(I0),
          R = floor(N0) - floor(S0) - floor(I0),
          Bcum = 0,
          Zcum = 0
        )
      })

      ## Transition events
      event_list <- list(
        c(S = 1, Bcum = 1),         # birth
        c(S = -1, I = 1, Zcum = 1), # infection
        c(I = -1, R = 1),           # removal
        c(S = -1),                  # natural mortality
        c(I = -1),
        c(R = -1)
      )

      ## Transition event rates
      compute_event_rates <- function(x, params, t) {
        m <- mu(t, params)
        c(
          nu(t, params),                     # birth
          beta(t, params) * x["S"] * x["I"], # infection
          gamma * x["I"] * (x["I"] > 1),     # removal
          m * x["S"],                        # natural mortality
          m * x["I"] * (x["I"] > 1),
          m * x["R"]
        )
      }

    } else if (model == "seir") {

      ## Initial state
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

      ## Transition events
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

      ## Transition event rates
      compute_event_rates <- function(x, params, t) {
        m <- mu(t, params)
        c(
          nu(t, params),                          # birth
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

    ## Generate a realization of the stochastic process
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

    ## For each time in `times`, find the index of the
    ## last previous time in `df$t`
    ind_state_out <- sapply(times, function(s) max(which(df$t <= s)))

    ## States at those times are states at times `times`
    df <- df[ind_state_out, ]
    df$t <- times
    rownames(df) <- NULL


  ### Simulate S(E)IR equations ... ------------------------------------
  ### if without demographic stochasticity

  } else {
    if (model == "sir") {

      ## Initial state
      x_init <- with(par_list, {
        c(
          S    = S0,
          logI = log(I0),
          R    = N0 - S0 - I0,
          Bcum = 0,
          Zcum = 0
        )
      })

      ## System of SIR equations
      compute_ode_rates <- function(t, y, parms) {
        m <- mu(t, parms)
        b <- beta(t, parms)
        y["I"] <- exp(y["logI"])
        dBcum <- nu(t, parms)
        dZcum <- b * y["S"] * y["I"]
        dS <- dBcum - dZcum - m * y["S"]
        dlogI <- b * y["S"] - gamma - m
        dR <- gamma * y["I"] - m * y["R"]
        list(c(dS, dlogI, dR, dBcum, dZcum))
      }

    } else if (model == "seir") {

      ## Initial state
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

      ## System of SEIR equations
      compute_ode_rates <- function(t, y, parms) {
        m <- mu(t, parms)
        b <- beta(t, parms)
        y["E"] <- exp(y["logE"])
        y["I"] <- exp(y["logI"])
        dBcum <- nu(t, parms)
        dZcum <- b * y["S"] * y["I"]
        dS <- dBcum - dZcum - m * y["S"]
        dlogE <- b * y["S"] * y["I"] / y["E"] - sigma - m
        dlogI <- sigma * y["E"] / y["I"] - gamma - m
        dR <- gamma * y["I"] - m * y["R"]
        list(c(dS, dlogE, dlogI, dR, dBcum, dZcum))
      }

    }

    ## Numerically integrate the system of S(E)IR equations
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


  ### Append other desired variables -----------------------------------

  df <- transform(df,
    t_years  = t * par_list$dt_days / 365,
    beta     = beta(t, par_list),
    mu       = mu(t, par_list),
    nu       = nu(t, par_list),
    p        = p,
    N        = S + (if (model == "seir") E else 0) + I + R,
    B        = c(NA, diff(Bcum)), # `B` from `Bcum` by first differences
    Z        = c(NA, diff(Zcum))  # `Z` from `Zcum` by first differences
  )
  df$B[1] <- df$B[2]
  df$Z[1] <- df$Z[2]


  ### Introduce observation error --------------------------------------

  ## `Zrep` from `Z` by binomial sampling
  b <- max(which(delay_dist > 0)) - 1
  Zrep <- rbinom(
    ## number of experiments
    n    = b + nrow(df),
    ## number of Bernoulli trials in each experiment
    size = round(c(rep(df$Z[1], b), df$Z)),
    ## success probability in each experiment
    prob = c(rep(df$p[1], b), df$p)
  )
  df$Zrep <- Zrep[b+(1:n)]

  ## `C` from `Zrep` by binning forward in time
  convol_out <- convol(Zrep, delay_dist, n = 1)
  df$C <- convol_out$simulation[b+(1:n), 1]

  if (model == "sir") {
    df <- df[, c("t", "t_years", "beta", "mu", "nu", "B",
                 "N", "S", "I", "Z", "p", "Zrep", "C")]
  } else if (model == "seir") {
    df <- df[, c("t", "t_years", "beta", "mu", "nu", "B",
                 "N", "S", "E", "I", "Z", "p", "Zrep", "C")]
  }
  structure(df, call = match.call(), arg_list = arg_list)
}
