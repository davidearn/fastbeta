make_data <- function(par_list       = list(),
                      n              = 20 * 365 / 7,
                      with_dem_stoch = TRUE,
                      seed           = NA,
                      ode_control    = list(
                        method = "lsoda",
                        rtol   = 1e-06,
                        atol   = 1e-06
                      )) {

## 1. Set-up -----------------------------------------------------------

# Create 3 `seeds` from 1 `seed`
if (!is.na(seed)) {
  set.seed(seed)
  seeds <- sample(1e09, 3)
}

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
t_out <- t0 + 0:n

# Gaussian white noise
if (!is.na(seed)) set.seed(seeds[1])
phi <- stats::rnorm(
  n    = length(t_out),
  mean = 0,
  sd   = epsilon
)

# Linear interpolant of `phi`
interpolate_phi <- stats::approxfun(
  x      = t_out,
  y      = phi,
  method = "linear",
  rule   = 2 # return `y[1]` and `y[length(y)]` outside range of `x`
)

# Seasonally forced transmission rate, without environmental noise
beta <- function(t) {
  beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
}

# Seasonally forced transmission rate, with environmental noise
beta_phi <- function(t) {
  beta_mean *
    (1 + alpha * cos(2 * pi * t / one_year + interpolate_phi(t)))
}


## 2.(a) Simulate SIR equations ... ------------------------------------
##       if with demographic stochasticity

if (with_dem_stoch) {

  ## NOTE: The adaptivetau package insists that simulations start at
  ##       time 0. To get simulations from time `t0`, we must take care
  ##       to add `t0` to the simulation time `t` where necessary.

  # Initial state
  x_init <- c(
    S = ceiling(S0),                           # susceptibles
    I = ceiling(I0),                           # infecteds
    R = round(N0) - ceiling(S0) - ceiling(I0), # removeds
    Q = 0                                      # cum. incidence
  )

  # Transition events
  event_list <- list(
    c(S = 1),                # birth
    c(S = -1, I = 1, Q = 1), # infection
    c(I = -1, R = 1),        # removal
    c(S = -1),               # natural mortality
    c(I = -1),
    c(R = -1)
  )

  # Transition event rates
  compute_event_rates <- function(x, params, t) {
    with(as.list(c(x, params)),
      {
        c(
          nu * hatN0,               # birth
          beta_phi(t + t0) * S * I, # infection
          gamma * I * (I > 1),      # removal
          mu * S,                   # natural mortality
          mu * I * (I > 1),
          mu * R
        )
      }
    )
  }

  # Generate a realization of the stochastic process
  if (!is.na(seed)) set.seed(seeds[2])
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
  colnames(df) <- c("t", "S", "I", "R", "Q")
  df <- transform(df, t = t + t0)

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
    S = ceiling(S0),                           # susceptibles
    logI = log(ceiling(I0)),                   # log infecteds
    R = round(N0) - ceiling(S0) - ceiling(I0), # removeds
    Q = 0                                      # cum. incidence
  )

  # System of SIRQ equations
  compute_sirq_rates <- function(t, y, parms) {
    with(as.list(c(y, parms)),
      {
        dS <- nu * hatN0 - beta_phi(t) * S * exp(logI) - mu * S
        dlogI <- beta_phi(t) * S - gamma - mu
        dR <- gamma * exp(logI) - mu * R
        dQ <- beta_phi(t) * S * exp(logI)
        list(c(dS, dlogI, dR, dQ))
      }
    )
  }

  # Create a list of arguments to be passed to `ode()`
  ode_args <- within(ode_control,
    {
      y     <- x_init
      times <- t_out
      func  <- compute_sirq_rates
      parms <- NULL # already in environment
    }
  )

  # Numerically integrate the system of SIRQ equations
  df <- as.data.frame(do.call(deSolve::ode, ode_args))
  colnames(df) <- c("t", "S", "logI", "R", "Q")
  df <- transform(df, I = exp(logI))

  # Warn if `ode()` returned early with unrecoverable error.
  # Append rows of `NA` until `nrow(df) = length(t_out)`.
  if (any(is.na(df))) {
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


## 3. Construct incidence `Z` and reported incidence `C` ... -----------
##    from cumulative incidence `Q` ------------------------------------

# `Z` from `Q` via first differences
df$Z <- c(NA, df$Q[-1] - df$Q[-nrow(df)])

# `C` from `Z` via delayed binomial sampling. `rbinom()`
# will warn about `NA` in the `size` argument. The warning
# can safely be suppressed.
reports_without_delay <- suppressWarnings(
  {
    if (!is.na(seed)) set.seed(seeds[3])
    stats::rbinom(
      n    = nrow(df),    # number of experiments
      size = round(df$Z), # number of Bernoulli trials
      p    = prep         # success probability
    )
  }
)
trepr <- round(trep)
reports_with_delay <- c(
  rep(NA, trepr),
  reports_without_delay[1:(nrow(df)-trepr)]
)
df$C <- reports_with_delay


## 4. Append other useful information ----------------------------------

df <- transform(df,
  t_years  = t * dt_weeks * (7 / 365),
  beta     = beta(t),
  beta_phi = beta_phi(t),
  N        = S + I + R
)


df <- df[, c("t", "t_years", "beta", "beta_phi",
             "N", "S", "I", "R", "Q", "Z", "C")]
attr(df, "arg_list") <-
  as.list(environment())[names(formals(make_data))]
df
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("dt_weeks", "t0", "prep", "trep", "hatN0", "N0", "S0", "I0",
      "nu", "mu", "tgen", "beta_mean", "alpha", "epsilon",
      "S", "logI", "R")
  )
}
