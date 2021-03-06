#' Simulate time series data
#'
#' `make_data()` simulates time series data according to a system of
#' SIR equations and a supplied list of parameter values. Observations
#' are recorded at equally spaced time points
#' \ifelse{latex}{\out{$t_k = t_0 + k \Delta t$}}{\ifelse{html}{\out{<i>t<sub>k</sub></i> = <i>t</i><sub>0</sub>+<i>k&Delta;t</i>}}{t_k = t_0 + k*Dt}}
#' (for \ifelse{latex}{\out{$k = 0,\ldots,n$}}{\ifelse{html}{\out{<i>k</i> = 0,...,<i>n</i>}}{k = 0,...,n}}),
#' where
#' \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#' denotes the observation interval.
#'
#' @section Model:
#' `make_data()` generates data at times
#' \ifelse{latex}{\out{$t_k = t_0 + k \Delta t$}}{\ifelse{html}{\out{<i>t<sub>k</sub></i> = <i>t</i><sub>0</sub>+<i>k&Delta;t</i>}}{t_k = t_0 + k*Dt}}
#' (for \ifelse{latex}{\out{$k = 0,\ldots,n$}}{\ifelse{html}{\out{<i>k</i> = 0,...,<i>n</i>}}{k = 0,...,n}})
#' using the system of SIR equations
#'
#' \ifelse{latex}{
#'   \out{
#'     \begin{array}[rlc]
#'       S' & = & \nu_\text{c} \hat{N}_0 - \beta(t) S I - \mu_\text{c} S \\
#'       I' & = & \beta(t) S I - \gamma I - \mu_\text{c} I \\
#'       R' & = & \gamma I - \mu_\text{c} R \\
#'     \end{array}
#'   }
#' }{
#'   \ifelse{html}{
#'     \out{
#'       <i>S</i>&prime; = <i>&nu;</i><sub>c</sub><i>&Ntilde;</i><sub>0</sub> &minus; <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&mu;</i><sub>c</sub><i>S</i><br>
#'       <i>I</i>&prime; = <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&gamma;I</i> &minus; <i>&mu;</i><sub>c</sub><i>I</i><br>
#'       <i>R</i>&prime; = <i>&gamma;I</i> &minus; <i>&mu;</i><sub>c</sub><i>R</i>
#'     }
#'   }{
#'     % S' = nu_c*hatN_0 - beta(t)*S*I - mu_c*S \cr
#'     % I' = beta(t)*S*I - gamma*I - mu_c*I \cr
#'     % R' = gamma*I - mu_c*R
#'   }
#' }
#'
#' with
#' \ifelse{latex}{\out{$\gamma = 1 / t_\text{gen}$}}{\ifelse{html}{\out{<i>&gamma;</i> = 1 / <i>t</i><sub>gen</sub>}}{gamma = 1/t_gen}}
#' and
#' \ifelse{latex}{\out{$\beta(t) = \langle\beta\rangle \left(1 + \alpha \cos\left(\frac{2 \pi t}{\text{1 year}}\right)\right)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>) = &langle;<i>&beta;</i>&rangle; (1 + <i>&alpha;</i> cos(2<i>&pi;t</i> / (1 year) + <i>&straightphi;</i>(<i>t</i>;<i>&straightepsilon;</i>)))}}{beta(t) = <beta>*(1 + alpha*cos(2*pi*t/(1 year)+phi(t;epsilon)))}},
#' where
#' \ifelse{latex}{\out{$\phi(t;\epsilon)$}}{\ifelse{html}{\out{<i>&straightphi;</i>(<i>t</i>;<i>&straightepsilon;</i>)}}{phi(t;epsilon)}}
#' is the linear interpolant of random noise
#' \ifelse{latex}{\out{$\{(t_k,\Phi_k)\}$}}{\ifelse{html}{\out{{(<i>t<sub>k</sub></i>,<i>&Phi;<sub>k</sub></i>)}}}{{(t_k,Phi_k)}}}
#' with
#' \ifelse{latex}{\out{$\Phi_k \sim \mathrm{Normal}(0,\epsilon^2)$}}{\ifelse{html}{\out{<i>&Phi;<sub>k</sub></i> ~ Normal(0,<i>&straightepsilon;</i><sup>2</sup>)}}{Phi_k ~ Normal(0,epsilon^2)}}.
#'
#' @section Randomness:
#' Simulations have three sources of randomness:
#'
#' 1. Environmental stochasticity
#'
#'    * The seasonal forcing function contains a randomly generated
#'      phase shift
#'      \ifelse{latex}{\out{$\phi(t;\epsilon)$}}{\ifelse{html}{\out{<i>&straightphi;</i>(<i>t</i>;<i>&straightepsilon;</i>)}}{phi(t;epsilon)}}.
#'      Randomness is introduced by choosing `epsilon` greater than 0.
#'
#' 2. Demographic stochasticity
#'
#'    * If `with_dem_stoch = TRUE`, then observations are generated
#'      by realizing a continuous-time stochastic process, in which
#'      event probabilities are determined by terms in the ODE.
#'      See [adaptivetau::ssa.adaptivetau()].
#'    * Probability of infected decrease is set to zero whenever the
#'      number of infecteds is 1, in order to prevent disease fadeout.
#'    * If `with_dem_stoch = FALSE`, then observations are generated
#'      by numerically integrating the ODE. See [deSolve::ode()].
#'
#' 3. Under-reporting of cases
#'
#'    * Reported incidence `C` is obtained from incidence `Z`
#'      by modeling `C[i]` as the number of successes in
#'      `Z[i-round(trep)]` independent Bernoulli trials,
#'      with success probability
#'      \ifelse{latex}{\out{$p_\text{rep}$}}{\ifelse{html}{\out{<i>p</i><sub>rep</sub>}}{p_rep}}.
#'      Randomness is introduced by choosing `prep` in (0,1).
#'
#' @param par_list A list containing:
#'
#'   \describe{
#'     \item{`dt_weeks`}{\[ \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}} \]
#'       Observation interval in weeks.
#'     }
#'     \item{`t0`}{\[ \ifelse{latex}{\out{$t_0$}}{\ifelse{html}{\out{<i>t</i><sub>0</sub>}}{t_0}} \]
#'       Time of the first observation in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     }
#'     \item{`prep`}{\[ \ifelse{latex}{\out{$p_\text{rep}$}}{\ifelse{html}{\out{<i>p</i><sub>rep</sub>}}{p_rep}} \]
#'       Case reporting probability.
#'     }
#'     \item{`trep`}{\[ \ifelse{latex}{\out{$t_\text{rep}$}}{\ifelse{html}{\out{<i>t</i><sub>rep</sub>}}{t_rep}} \]
#'       Case reporting delay in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     }
#'     \item{`hatN0`}{\[ \ifelse{latex}{\out{$\widehat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}} \]
#'       Population size at time
#'       \ifelse{latex}{\out{$t = 0$}}{\ifelse{html}{\out{<i>t</i> = 0}}{t = 0}}
#'       years.
#'     }
#'     \item{`N0`}{\[ \ifelse{latex}{\out{$N_0$}}{\ifelse{html}{\out{<i>N</i><sub>0</sub>}}{N_0}} \]
#'       Population size at time
#'       \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'     }
#'     \item{`S0`}{\[ \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}} \]
#'       Number of susceptibles at time
#'       \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'     }
#'     \item{`I0`}{\[ \ifelse{latex}{\out{$I_0$}}{\ifelse{html}{\out{<i>I</i><sub>0</sub>}}{I_0}} \]
#'       Number of infecteds at time
#'       \ifelse{latex}{\out{$t = t_0$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub>0</sub>}}{t = t_0}}.
#'     }
#'     \item{`nu`}{\[ \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}} \]
#'       Birth rate expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       and relative to
#'       \ifelse{latex}{\out{$\hat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}}.
#'     }
#'     \item{`mu`}{\[ \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}} \]
#'       Natural mortality rate expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       and per capita.
#'     }
#'     \item{`tgen`}{\[ \ifelse{latex}{\out{$t_\text{gen}$}}{\ifelse{html}{\out{<i>t</i><sub>gen</sub>}}{t_gen}} \]
#'       Mean generation interval of the disease of interest in units
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     }
#'     \item{`beta_mean`}{\[ \ifelse{latex}{\out{$\langle\beta\rangle$}}{\ifelse{html}{\out{&langle;<i>&beta;</i>&rangle;}}{<beta>}} \]
#'       Mean of the seasonally forced transmission rate
#'       \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#'       expressed per unit
#'       \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'       per susceptible per infected.
#'     }
#'     \item{`alpha`}{\[ \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}} \]
#'       Amplitude of the seasonally forced transmission rate
#'       \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#'       relative to the mean.
#'     }
#'     \item{`epsilon`}{\[ \ifelse{latex}{\out{$\epsilon$}}{\ifelse{html}{\out{<i>&straightepsilon;</i>}}{epsilon}} \]
#'       Standard deviation of the random phase shift in the seasonally
#'       forced transmission rate
#'       \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}.
#'     }
#'   }
#' @param n Numeric scalar. Time between the first and last observations
#'   in units
#'   \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}},
#'   so that simulated time series have `n+1` observations.
#' @param with_dem_stoch Logical scalar. If `TRUE`, then the simulation
#'   is generated by [adaptivetau::ssa.adaptivetau()]. Otherwise, it is
#'   generated by [deSolve::ode()] (see Details).
#' @param seed Integer scalar or `NA`. If an integer, then `seed`
#'   is passed to `set.seed()` internally, making the simulation,
#'   which is randomly generated (see Details), reproducible with
#'   an identical call to `make_data()`.
#' @param ode_control A list of arguments to be passed to
#'   [deSolve::ode()], specifying options for numerical integration,
#'   such as `method`, `rtol`, and `atol`. Not used if
#'   `with_dem_stoch = TRUE`.
#'
#' @return
#' A data frame with numeric columns:
#'
#' \describe{
#'   \item{`t`}{Time in units
#'     \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}.
#'     Equal to `t0 + seq(0, n, by = 1)`.
#'   }
#'   \item{`t_years`}{Time in years. Equal to
#'     `(t0 + seq(0, n, by = 1)) * dt_weeks * (7 / 365)`.
#'   }
#'   \item{`beta`}{Seasonally forced transmission rate expressed per
#'     unit
#'     \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}},
#'     per susceptible per infected, without environmental noise.
#'   }
#'   \item{`beta_phi`}{Seasonally forced transmission rate expressed per
#'     unit
#'     \ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
#'     per susceptible per infected, with environmental noise.
#'   }
#'   \item{`N`}{Population size.}
#'   \item{`S`}{Number of susceptibles.}
#'   \item{`I`}{Number of infecteds.}
#'   \item{`R`}{Number of removeds.}
#'   \item{`Q`}{Cumulative incidence. `Q[i]` is the number of infections
#'     between times `t[1]` and `t[i]`.
#'   }
#'   \item{`Z`}{Incidence. `Z[i]` is the number of infections
#'     between times `t[i-1]` and `t[i]`, computed as `Q[i] - Q[i-1]`.
#'   }
#'   \item{`C`}{Reported incidence. `C[i]` is the number of infections
#'     reported between times `t[i-1]` and `t[i]`, equal to the number
#'     of successes in `Z[i-round(trep)]` independent Bernoulli trials,
#'     with success probability
#'     \ifelse{latex}{\out{$p_\text{rep}$}}{\ifelse{html}{\out{<i>p</i><sub>rep</sub>}}{p_rep}}.
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
#'   n = 20 * 365 / 7,       # 20 years is ~1042 weeks
#'   with_dem_stoch = FALSE  # numerical integration of d(S,I,R)/dt
#' )
#' head(df)
#' 
#' # Reproducible stochastic simulation
#' par_list <- make_par_list(
#'   dt_weeks = 1,
#'   epsilon = 0.5, # s.d. of environmental noise
#'   prep = 0.5     # case reporting probability
#' )
#' df <- make_data(
#'   par_list = par_list,
#'   n = 20 * 365 / 7,      # 20 years is ~1042 weeks
#'   with_dem_stoch = TRUE, # stochastic simulation of d(S,I,R)/dt
#'   seed = 5               # seed for RNG
#' )
#' head(df)
#'
#' @md
#' @export
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
