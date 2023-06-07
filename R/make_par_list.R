make_par_list <- function(dt_weeks  = 1,
                          t0        = 2000 * (365 / 7) / dt_weeks,
                          prep      = 1,
                          trep      = 0,
                          hatN0     = 1e06,
                          N0        = NA,
                          S0        = NA,
                          I0        = NA,
                          nu        = 0.04 * (7 / 365) * dt_weeks,
                          mu        = 0.04 * (7 / 365) * dt_weeks,
                          tgen      = 13 * (1 / 7) / dt_weeks,
                          Rnaught   = 20,
                          beta_mean = NA,
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

# If `Rnaught` was defined
if (!is.na(Rnaught)) {
  beta_mean <- (mu / (nu * hatN0)) * Rnaught * (gamma + mu)
# If `Rnaught` was not defined but `beta_mean` was
} else if (is.na(Rnaught) && !is.na(beta_mean)) {
  Rnaught <- ((nu * hatN0) / mu) * beta_mean / (gamma + mu)
# If neither was defined
} else {
  stop(
    "At most one of `Rnaught` and `beta_mean` can be `NA`.",
    call. = FALSE
  )
}


# If `N0`, `S0`, or `I0` is not defined, then
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

  # Seasonally forced transmission rate
  beta <- function(t) {
    beta_mean * (1 + alpha * cos(2 * pi * t / one_year))
  }

  # System of SIR equations
  compute_sir_rates <- function(t, y, parms) {
    with(as.list(c(y, parms)),
      {
        dS <- nu * hatN0 - beta(t) * S * exp(logI) - mu * S
        dlogI <- beta(t) * S - gamma - mu
        dR <- gamma * exp(logI) - mu * R
        list(c(dS, dlogI, dR))
      }
    )
  }

  # List of arguments to be passed to `ode()`
  ode_args <- within(ode_control,
    {
      y     <- x_init
      times <- t_out
      func  <- compute_sir_rates
      parms <- NULL # already in environment
    }
  )

  # Numerically integrate the system of SIR equations
  df <- as.data.frame(do.call(deSolve::ode, ode_args))

  # Assign final values of `S+I+R`, `S`, and `I`
  if (is.na(N0)) {
    N0 <- sum(df[nrow(df), c("S", "R")]) + exp(df[nrow(df), "logI"])
  }
  if (is.na(S0)) {
    S0 <- df[nrow(df), "S"]
  }
  if (is.na(I0)) {
    I0 <- exp(df[nrow(df), "logI"])
  }

  # Warn if `ode()` returned early with unrecoverable error
  if (any(is.na(df))) {
    warning(
      "`ode()` could not complete the integration. ",
      "Retry with modified `ode_control`.",
      call. = FALSE
    )
  }
}


out <- as.list(environment())[names(formals(make_par_list))]
out$ode_control <- NULL
out
}
