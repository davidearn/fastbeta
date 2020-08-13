#' \loadmathjax
#' Simulate epidemic time series data
#'
#' @description
#' `make_data()` simulates epidemic time series data using a system of
#' SIR equations and a supplied list of parameter values. Observations
#' are recorded at equally spaced time points \mjseqn{t_i = t_0 + i \Delta t}
#' (for \mjseqn{i = 0, \ldots, n}). Among other things, the simulation model
#' assumes:
#' * A seasonally forced transmission rate.
#' * Constant birth and per capita natural mortality rates.
#' * Infections reported with a fixed probability after a
#'   delay that either (i) is fixed or (ii) follows a
#'   negative binomial distribution.
#'
#' @param par_list A list of parameter values with elements:
#'   \describe{
#'     \item{`dt_weeks`}{\mjseqn{\lbrace\,\Delta t\,\rbrace}
#'       Numeric scalar. Observation interval in weeks.
#'     }
#'     \item{`t0`}{\mjseqn{\lbrace\,t_0/\Delta t\,\rbrace}
#'       Numeric scalar. Time of the first observation
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`prep`}{\mjseqn{\lbrace\,p_\text{rep}\,\rbrace}
#'       Numeric scalar. Probability that an infection is reported.
#'     }
#'     \item{`trep`}{\mjseqn{\lbrace\,t_\text{rep}/\Delta t\,\rbrace}
#'       Numeric scalar. Mean time between infection and reporting
#'       in units \mjseqn{\Delta t}.
#'     }
#'     \item{`k`}{\mjseqn{\lbrace\,k\,\rbrace}
#'       Numeric scalar. An optional dispersion parameter.
#'       If specified, then the time between infection and reporting
#'       in units \mjseqn{\Delta t} will be modeled as a negative
#'       binomially distributed random variable with dispersion `k`
#'       and mean `trep`. If not, then the delay will be fixed equal
#'       to `round(trep)`.
#'     }
#'     \item{`hatN0`}{\mjseqn{\lbrace\,\out{\widehat{N}_0}\,\rbrace}
#'       Numeric scalar. Population size at time \mjseqn{t = 0} years.
#'     }
#'     \item{`N0`}{\mjseqn{\lbrace\,N_0\,\rbrace}
#'       Numeric scalar. Population size at time \mjseqn{t = t_0}.
#'     }
#'     \item{`S0`}{\mjseqn{\lbrace\,S_0\,\rbrace}
#'       Numeric scalar. Number of susceptibles at time \mjseqn{t = t_0}.
#'     }
#'     \item{`I0`}{\mjseqn{\lbrace\,I_0\,\rbrace}
#'       Numeric scalar. Number of infecteds at time \mjseqn{t = t_0}.
#'     }
#'     \item{`nu`}{\mjseqn{\lbrace\,\nu \Delta t\,\rbrace}
#'       Numeric scalar. Birth rate expressed
#'       per unit \mjseqn{\Delta t} and relative to
#'       \mjseqn{\out{\widehat{N}_0}}.
#'     }
#'     \item{`mu`}{\mjseqn{\lbrace\,\mu \Delta t\,\rbrace}
#'       Numeric scalar. Natural mortality rate expressed
#'       per unit \mjseqn{\Delta t} and per capita.
#'     }
#'     \item{`tgen`}{\mjseqn{\lbrace\,t_\text{gen}/\Delta t\,\rbrace}
#'       Numeric scalar. Mean generation interval of the disease
#'       of interest in units \mjseqn{\Delta t}.
#'     }
#'     \item{`beta_mean`}{\mjseqn{\lbrace\,\langle\beta\rangle \Delta t\,\rbrace}
#'       Numeric scalar. Mean (long-term average) of the seasonally
#'       forced transmission rate \mjseqn{\beta(t)} expressed per unit
#'       \mjseqn{\Delta t} per susceptible per infected.
#'     }
#'     \item{`alpha`}{\mjseqn{\lbrace\,\alpha\,\rbrace}
#'       Numeric scalar. Amplitude of the seasonally forced
#'       transmission rate \mjseqn{\beta(t)} relative to the mean.
#'     }
#'     \item{`epsilon2`}{\mjseqn{\lbrace\,\epsilon^2\,\rbrace}
#'       Numeric scalar. Variance of the normally distributed
#'       phase shift introduced to the seasonally forced
#'       transmission rate \mjseqn{\beta(t)}.
#'     }
#'   }
#' @param n Integer scalar. Time between the first and last observations
#'   in units \mjseqn{\Delta t}, so that simulated time series have `n + 1`
#'   observations. (If numeric but not integer, then `n` is replaced by
#'   `floor(n)`.)
#' @param with_dem_stoch Logical scalar. If `TRUE`, then the simulation
#'   is generated using [adaptivetau::ssa.adaptivetau()]. Otherwise, it
#'   is generated using [deSolve::ode()] (see Details).
#'
#' @return
#' A data frame with `n + 1` rows corresponding to equally spaced time
#' points \mjseqn{t_i = t_0 + i \Delta t} (for \mjseqn{i = 0, \ldots, n}),
#' and numeric columns:
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
#'     per unit \mjseqn{\Delta t} per susceptible per infected,
#'     with environmental noise.
#'   }
#'   \item{`N`}{Population size.}
#'   \item{`S`}{Number of susceptibles.}
#'   \item{`I`}{Number of infecteds.}
#'   \item{`Z`}{Incidence. `Z[i]` is the number of infections
#'     between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`Zrep`}{Incidence conditional on reporting.
#'     `Zrep[i]` is the number of infections between times `t[i-1]` and `t[i]`
#'     that are eventually reported.
#'   }
#'   \item{`C`}{Reported incidence. `C[i]` is the number of infections
#'     reported between times `t[i-1]` and `t[i]`.
#'   }
#'   \item{`B`}{Births. `B[i]` is the number of births between times
#'     `t[i-1]` and `t[i]`. Equal to `with(par_list, rep(nu * hatN0, n + 1))`
#'     if `with_dem_stoch = FALSE`.
#'   }
#'   \item{`mu`}{Per capita natural mortality rate. `mu[i]` is rate
#'     at time `t[i]`. Equal to `with(par_list, rep(mu, n + 1))`.
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
#' and \mjseqn{\phi(t;\epsilon^2)} is the linear interpolant of noise
#' \mjseqn{\lbrace(t_k;\Phi_k)\rbrace_{k=0}^n} with
#' 
#' \mjsdeqn{\Phi_k \sim \mathrm{Normal}(0,\epsilon^2)\,,}
#'
#' modeling **environmental stochasticity**.
#'
#' `make_data()` generates observations of the above system at equally spaced
#' times \mjseqn{t_i = t_0 + i \Delta t} (for \mjseqn{i = 0, \ldots, n}) by
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
#' \mjseqn{Z_\text{cum}} via first differences:
#'
#' \mjsdeqn{\begin{align} B(t_i) &= B_\text{cum}(t_i) - B_\text{cum}(t_{i-1}), \cr Z(t_i) &= Z_\text{cum}(t_i) - Z_\text{cum}(t_{i-1}). \end{align}}
#'
#' It then does binomial sampling to generate, for each observation
#' of \mjseqn{Z}, a number of infections that is eventually reported:
#'
#' \mjsdeqn{Z_\text{rep}(t_i) \sim \mathrm{Binomial}\big(Z(t_i),p_\text{rep}\big)\,.}
#' 
#' Finally, it generates reported incidence \mjseqn{C} by sampling
#' from a reporting delay distribution and binning infections counted
#' by \mjseqn{Z_\text{rep}} forward in time. Compared to \mjseqn{Z},
#' \mjseqn{C} models **observation error**, i.e., under-reporting of
#' cases with a delay between infection and reporting. If `par_list$k`
#' is non-`NULL`, then the reporting delay in units \mjseqn{\Delta t}
#' is modeled as
#'
#' \mjsdeqn{D \sim \mathrm{NegativeBinomial}\left(k,\frac{k}{m+k}\right)\,,}
#'
#' where \mjseqn{k} is a dispersion parameter and
#' \mjseqn{m = t_\text{rep}/\Delta t} is the mean reporting delay.
#' In this case, infections counted in \mjseqn{Z_\text{rep}(t_i)}
#' are reported in \mjseqn{C(t_{i+D})}. On the other hand, if
#' `par_list$k` is `NULL`, then the reporting delay is fixed equal
#' to \mjseqn{d = \mathrm{nint}(t_\text{rep}/\Delta t)}, and so
#' \mjseqn{C(t_{i+d}) = Z_\text{rep}(t_i)\,.}
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
#'   epsilon2 = 1, # environmental sxtochasticity
#'   prep = 0.5, # random under-reporting of infections
#'   trep = 2, # random reporting delays
#'   k = 4
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
                      with_dem_stoch = TRUE) {
  ## 1. Set-up -----------------------------------------------------------

  # Save arguments in a list
  arg_list <- as.list(environment())

  # Load necessary elements of `par_list` into the execution environment
  list2env(
    par_list[c("dt_weeks", "t0", "prep", "trep", "k",
               "hatN0", "N0", "S0", "I0", "nu", "mu", "tgen",
               "beta_mean", "alpha", "epsilon2")],
    envir = environment()
  )

  # Some derived quantities
  gamma <- 1 / tgen
  one_year <- (365 / 7) / dt_weeks

  # Time points
  n <- floor(n)
  t_out <- t0 + 0:n

  # Environmental noise
  phi <- stats::rnorm(
    n    = length(t_out),
    mean = 0,
    sd   = sqrt(epsilon2)
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
    df <- as.data.frame(
      adaptivetau::ssa.adaptivetau(
        x_init, event_list, compute_event_rates,
        params    = NULL, # in enclosing environment of `compute_event_rates()`
        tf        = n,    # final time point
        tl.params = list( # other instructions:
          epsilon     = 0.05,
          delta       = 0.05,
          maxtau      = 0.5, # adaptive time step must not exceed 1
          extraChecks = TRUE
        )
      )
    )
    colnames(df)[1] <- "t"
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

    # Numerically integrate the system of SIR equations
    df <- as.data.frame(
      deSolve::ode(
        y     = x_init,
        times = t_out,
        func  = compute_sir_rates,
        parms = NULL # in enclosing environment of `compute_sir_rates()`
      )
    )
    colnames(df)[1] <- "t"
    df$I <- exp(df$logI)

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


  ## 4. Introduce observation error --------------------------------------

  # `Zrep` from `Z` via binomial sampling
  df$Zrep <- c(NA,
    stats::rbinom(
      n    = nrow(df) - 1,    # number of experiments
      size = round(df$Z[-1]), # numbers of Bernoulli trials
      p    = prep             # success probability
    )
  )

  # `C` from `Zrep` with a fixed shift forward
  if (!exists("k", inherits = FALSE) || is.null(k)) {
    trepr <- round(trep)
    df$C <- c(rep(NA, trepr), df$Zrep[1:(nrow(df)-trepr)])
  # `C` from `Zrep` via negative binomial sampling
  } else {
    counts <- df$Zrep[-1]
    delays <- stats::rnbinom(sum(counts), size = k, mu = trep)
    index_counts <- rep(seq_along(counts), times = counts)
    index_counts_with_delays <- index_counts + delays
    counts_with_delays <- rep(0, max(index_counts_with_delays))
    counts_with_delays[sort(unique(index_counts_with_delays))] <-
      table(index_counts_with_delays)
    df$C <- c(NA, counts_with_delays[1:(nrow(df)-1)])
  }
    
  df <- df[, c("t", "t_years", "beta", "beta_phi",
               "N", "S", "I", "Z", "Zrep", "C", "B", "mu")]
  attr(df, "call") <- match.call()
  attr(df, "arg_list") <- arg_list
  df
}

# For `R CMD check`
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("dt_weeks", "t0", "prep", "trep", "k",
      "hatN0", "N0", "S0", "I0", "nu", "mu", "tgen",
      "beta_mean", "alpha", "epsilon2",
      "S", "R", "Zcum", "Bcum")
  )
}
