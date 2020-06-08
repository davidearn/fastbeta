#' \loadmathjax
#' Simulate epidemic time series data
#'
#' `make_data()` simulates epidemic time series data using a system of
#' SIR equations and a supplied list of parameter values. Observations
#' are recorded at equally spaced time points
#' \mjeqn{t_k = t_0 + k \Delta t}{t_k = t_0 + k dt}
#' (for \mjeqn{k = 0, \ldots, n}{k = 0,...,n}).
#' Among other things, the simulation model assumes:
#' * A seasonally forced transmission rate.
#' * Constant birth and per capita natural mortality rates. 
#' * Cases reported after a fixed delay since infection,
#'   with a fixed probability.
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
#' \mjdeqn{\out{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu \widehat{N}_0 - \beta_\phi(t) S I - \mu S, \cr \frac{\text{d}I}{\text{d}t} &= \beta_\phi(t) S I - \gamma I - \mu I, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu R, \cr \frac{\text{d}B_\text{cum}}{\text{d}t} &= \nu \widehat{N}_0, \cr \frac{\text{d}Z_\text{cum}}{\text{d}t} &= \beta_\phi(t) S I. \end{align}}}{(1) dS/dt = nu hatN_0 - beta_phi(t) S I - mu S,    (2) dI/dt = beta_phi(t) S I - gamma I - mu I,    (3) dR/dt = gamma I - mu R,    (4) dB_cum/dt = nu hatN_0,    (5) dZ_cum/dt = beta_phi(t).}
#' 
#' Here, \mjeqn{\gamma = 1 / t_\text{gen}}{gamma = 1 / t_gen},
#' 
#' \mjdeqn{\beta_\phi(t) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}} + \phi(t;\epsilon)\right) \right\rbrack,}{beta(t) = <beta> (1 + alpha cos(2 pi t / (1 year) + phi(t;epsilon))),}
#' 
#' and \mjeqn{\phi(t;\epsilon)}{phi(t;epsilon)}
#' is the linear interpolant of noise
#' \mjeqn{\lbrace(t_k;\Phi_k)\rbrace_{k=0}^n}{{(t_k,Phi_k)}}
#' with
#' 
#' \mjdeqn{\Phi_k \sim \mathrm{Normal}(0,\epsilon^2),}{Phi_k ~ Normal(0,epsilon^2),}
#'
#' modeling **environmental stochasticity**.
#'
#' `make_data()` generates observations of the above system at equally
#' spaced times \mjeqn{t_k = t_0 + k \Delta t}{t_k = t_0 + k dt} (for
#' \mjeqn{k = 0, \ldots, n}{k = 0,...,n}) by either (i) numerically
#' integrating the ODE using [deSolve::ode()] (`with_dem_stoch = FALSE`),
#' or (ii) realizing a corresponding continuous-time stochastic process
#' using [adaptivetau::ssa.adaptivetau()] (`with_dem_stoch = TRUE`).
#' Both methods use initial state
#'
#' \mjdeqn{\begin{bmatrix} S(t_0) \cr I(t_0) \cr R(t_0) \cr B_\text{cum}(t_0) \cr Z_\text{cum}(t_0) \end{bmatrix} = \begin{bmatrix} S_0 \cr I_0 \cr N_0 - S_0 - I_0 \cr 0 \cr 0 \end{bmatrix}.}{(S(t_0),I(t_0),R(t_0),B_cum(t_0),Z_cum(t_0)) = (S_0,I_0,N_0 - S_0 - I_0,0,0).}
#' 
#' The latter method defines event probabilities as proportional to
#' terms in the ODE and models **demographic stochasticity**. Disease
#' fadeout in simulations with demographic stochasticity is prevented
#' by setting the probability of infected decrease to zero whenever
#' the number of infecteds is 1.
#'
#' `make_data()` calculates births \mjseqn{B} and incidence \mjseqn{Z}
#' from cumulative births \mjeqn{B_\text{cum}}{B_cum} and cumulative
#' incidence \mjeqn{Z_\text{cum}}{Z_cum} via first differences,
#'
#' \mjdeqn{\begin{align} B(t_k) &= B_\text{cum}(t_k) - B_\text{cum}(t_{k-1}), \cr Z(t_k) &= Z_\text{cum}(t_k) - Z_\text{cum}(t_{k-1}), \end{align}}{(1) B(t_k) = B_cum(t_k) - B_cum(t_{k-1}),    (2) Z(t_k) = Z_cum(t_k) - Z_cum(t_{k-1}),}
#'
#' then simulates reported incidence \mjseqn{C} from
#' incidence \mjseqn{Z} via lagged binomial sampling,
#'
#' \mjdeqn{C(t_{k+r}) \sim \mathrm{Binomial}\big(Z(t_k),p_\text{rep}\big),}{C(t_{k+r}) ~ Binomial(Z(t_k), p_rep),}
#'
#' modeling **observation error** (random under-reporting of cases with
#' a delay between infection and reporting). Here,
#' \mjeqn{r = \mathrm{nint}(t_\text{rep} / \Delta t)}{r = nint(t_rep / dt)}.
#'
#' @param par_list A list of parameter values containing:
#'   \describe{
#'     \item{`dt_weeks`}{\mjeqn{\lbrack \Delta t \rbrack}{\[dt\]}
#'       Numeric scalar. Observation interval in weeks.
#'     }
#'     \item{`t0`}{\mjeqn{\lbrack t_0 \rbrack}{\[t_0\]}
#'       Numeric scalar. Time of the first observation
#'       in units \mjeqn{\Delta t}{dt}.
#'     }
#'     \item{`prep`}{\mjeqn{\lbrack p_\text{rep} \rbrack}{\[p_rep\]}
#'       Numeric scalar. Case reporting probability.
#'     }
#'     \item{`trep`}{\mjeqn{\lbrack t_\text{rep} \rbrack}{\[t_rep\]}
#'       Numeric scalar. Case reporting delay in units
#'       \mjeqn{\Delta t}{dt}.
#'     }
#'     \item{`hatN0`}{\mjeqn{\lbrack \out{\widehat{N}_0} \rbrack}{\[hatN_0\]}
#'       Numeric scalar. Population size at time \mjseqn{t = 0} years.
#'     }
#'     \item{`N0`}{\mjeqn{\lbrack N_0 \rbrack}{\[N_0\]}
#'       Numeric scalar. Population size at time \mjseqn{t = t_0}.
#'     }
#'     \item{`S0`}{\mjeqn{\lbrack S_0 \rbrack}{\[S_0\]}
#'       Numeric scalar. Number of susceptibles at time \mjseqn{t = t_0}.
#'     }
#'     \item{`I0`}{\mjeqn{\lbrack I_0 \rbrack}{\[I_0\]}
#'       Numeric scalar. Number of infecteds at time \mjseqn{t = t_0}.
#'     }
#'     \item{`nu`}{\mjeqn{\lbrack \nu_\text{c} \rbrack}{\[nu_c\]}
#'       Numeric scalar. Birth rate expressed
#'       per unit \mjeqn{\Delta t}{dt} and
#'       relative to \mjeqn{\out{\widehat{N}_0}}{hatN_0}.
#'     }
#'     \item{`mu`}{\mjeqn{\lbrack \mu_\text{c} \rbrack}{\[mu_c\]}
#'       Numeric scalar. Natural mortality rate expressed
#'       per unit \mjeqn{\Delta t}{dt} and per capita.
#'     }
#'     \item{`tgen`}{\mjeqn{\lbrack t_\text{gen} \rbrack}{\[t_gen\]}
#'       Numeric scalar. Mean generation interval of the disease
#'       of interest in units \mjeqn{\Delta t}{dt}.
#'     }
#'     \item{`beta_mean`}{\mjeqn{\lbrack \langle\beta\rangle \rbrack}{\[<beta>\]}
#'       Numeric scalar. Mean (long-term average) of the seasonally
#'       forced transmission rate \mjeqn{\beta(t)}{beta(t)} expressed
#'       per unit \mjeqn{\Delta t}{dt} per susceptible per infected.
#'     }
#'     \item{`alpha`}{\mjeqn{\lbrack \alpha \rbrack}{\[alpha\]}
#'       Numeric scalar. Amplitude of the seasonally forced transmission
#'       rate \mjeqn{\beta(t)}{beta(t)} relative to the mean.
#'     }
#'     \item{`epsilon`}{\mjeqn{\lbrack \epsilon \rbrack}{\[epsilon\]}
#'       Numeric scalar. Standard deviation of the random phase shift
#'       introduced to the seasonally forced transmission rate
#'       \mjeqn{\beta(t)}{beta(t)}.
#'     }
#'   }
#' @param n Integer scalar. Time between the first and last observations
#'   in units \mjeqn{\Delta t}{dt}, so that simulated time series have
#'   `n + 1` observations. (If numeric but not integer, then `n` is
#'   replaced by `floor(n)`.)
#' @param with_dem_stoch Logical scalar. If `TRUE`, then the simulation
#'   is generated using [adaptivetau::ssa.adaptivetau()]. Otherwise, it
#'   is generated using [deSolve::ode()] (see Details).
#' @param ode_control A list of optional arguments of [deSolve::ode()],
#'   specifying options for numerical integration, such as `method`,
#'   `rtol`, and `atol`. Not used if `with_dem_stoch = TRUE` (see Details).
#'
#' @return
#' A data frame with `n + 1` rows corresponding to equally spaced times
#' \mjeqn{t_k = t_0 + k \Delta t}{t_k = t_0 + k dt}
#' (for \mjeqn{k = 0, \ldots, n}{k = 0,...,n}), and numeric columns:
#'
#' \describe{
#'   \item{`t`}{Time in units \mjeqn{\Delta t}{dt}. Equal to
#'     `par_list$t0 + 0:n`.
#'   }
#'   \item{`t_years`}{Time in years. Equal to
#'     `(par_list$t0 + 0:n) * par_list$dt_weeks * (7 / 365)`.
#'   }
#'   \item{`beta`}{Seasonally forced transmission rate expressed
#'     per unit \mjeqn{\Delta t}{dt} per susceptible per infected,
#'     without environmental noise.
#'   }
#'   \item{`beta_phi`}{Seasonally forced transmission rate expressed
#'     per unit \mjeqn{\Delta t}{dt} per susceptible per infected,
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
#'     for convenience, mainly because [estimate_beta_SI()] requires
#'     a data frame with columns `C`, `B`, and `mu`.
#'   }
#' }
#'
#' A list of the arguments of `make_data()` is included as an attribute.
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
#' @md
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
      "`ode()` could not complete the integration. ",
      "Retry with modified `ode_control`.",
      call. = FALSE
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
attr(df, "arg_list") <-
  as.list(environment())[names(formals(make_data))]
df
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("dt_weeks", "t0", "prep", "trep", "hatN0", "N0", "S0", "I0",
      "nu", "mu", "tgen", "beta_mean", "alpha", "epsilon",
      "S", "R", "Zcum", "Bcum")
  )
}
