estimate_beta_SI <- function(df       = data.frame(),
                             par_list = list(),
                             method   = c("trapezoid", "both")) {

## 1. Set-up -----------------------------------------------------------

# Load necessary elements of `par_list` into the execution environment
list2env(
  par_list[c("prep", "trep", "S0", "I0", "hatN0", "nu", "mu", "tgen")],
  envir = environment()
)

# Logical: Is `df` missing a column `B`? What about `mu`?
B_was_not_def <- is.null(df$B)
mu_was_not_def <- is.null(df$mu)

# Preallocate memory for output. Assume constant vital rates,
# if necessary.
df <- data.frame(
  t     = df$t,
  C     = df$C,
  Z     = NA,
  B     = if (B_was_not_def) nu * hatN0 * 1 else df$B,
  mu    = if (mu_was_not_def) mu else df$mu,
  S     = NA,
  I     = NA,
  beta  = NA
)


## 2. Impute missing values in reported incidence, ... -----------------
##    births, and natural mortality rate

df[c("C", "B", "mu")] <- lapply(df[c("C", "B", "mu")],
  function(x) {
    if (!any(is.na(x))) {
      return(x)
    }

    # Indices of missing values
    ind_na <- which(is.na(x))

    # Function that performs linear interpolation
    # between time points where `x` is observed
    impute_na <- stats::approxfun(
      x      = df$t[-ind_na],
      y      = x[-ind_na],
      method = "linear",
      rule   = 1 # return `NA` outside range of (argument) `x`
    )

    # Impute missing values. Missing values outside of
    # the range of observed data are retained as `NA`.
    replace(x, ind_na, impute_na(df$t[ind_na]))
  }
)


## 3. Impute zeros in reported incidence -------------------------------

df$C <- local(
  {
    if (!any(df$C == 0, na.rm = TRUE)) {
      return(df$C)
    }

    # Indices of positive values
    ind_pos <- which(df$C > 0)

    # Indices of zeros between two positive values
    ind_zero <- which(df$C == 0)
    ind_zero <- ind_zero[
      ind_zero > min(ind_pos) &
      ind_zero < max(ind_pos)
    ]

    # Function that performs linear interpolation
    # between time points where `df$C` is positive
    impute_zero <- stats::approxfun(
      x      = df$t[ind_pos],
      y      = df$C[ind_pos],
      method = "linear",
      rule   = 1 # return `NA` outside range of `x`
    )

    # Impute zeros. Zeros outside of the range of
    # positive data are retained as zeros due to
    # the definition of `ind_zero`.
    replace(df$C, ind_zero, impute_zero(df$t[ind_zero]))
  }
)


## NOTE: In the best case, there are no longer missing values or zeros
##       in `df$C`. In the worst case, `df$C` now looks like
##       `c(NA,...,NA,0,...,0,+,...,+,0,...,0,NA,...,NA)`,
##       where `+` denotes a positive number.


## 4. Estimate incidence -----------------------------------------------

trepr <- round(trep)
df$Z <- c(
  (1 / prep) * df$C[(trepr+1):nrow(df)],
  rep(NA, trepr)
)


## 5. Estimate susceptibles and infecteds ------------------------------

gamma <- 1 / tgen

df[c("S", "I")] <- with(df,
  {
    S[1] <- S0
    I[1] <- I0
    if (method[1] == "trapezoid") {
      for (i in 2:nrow(df)) {
        S[i] <- ((1 - 0.5 * mu[i-1] * 1) * S[i-1] + B[i] - Z[i]) /
          (1 + 0.5 * mu[i] * 1)
        I[i] <- ((1 - 0.5 * (gamma + mu[i-1]) * 1) * I[i-1] + Z[i]) /
          (1 + 0.5 * (gamma + mu[i]) * 1)
      }
      beta <- (Z + c(Z[-1], NA)) / (2 * S * I * 1)
    } else if (method[1] == "backward") {
      for (i in 2:nrow(df)) {
        S[i] <- (S[i-1] + B[i] - Z[i]) / (1 + mu[i] * 1)
        I[i] <- (I[i-1] + Z[i]) / (1 + (gamma + mu[i]) * 1)
      }
      beta <- Z / (S * I * 1)
    } else if (method[1] == "forward") {
      for (i in 2:nrow(df)) {
        S[i] <- (1 - mu[i-1] * 1) * S[i-1] + B[i] - Z[i]
        I[i] <- (1 - (gamma - mu[i-1]) * 1) * I[i-1] + Z[i]
      }
      beta <- c(Z[-1], NA) / (S * I * 1)
    } else {
      stop(
        "SI method: `method[1]` must be ",
        "\"forward\", \"backward\", or \"trapezoid\".",
        call. = FALSE
      )
    }
    list(S, I)
  }
)


## 6. Estimate transmission rate ---------------------------------------

df$beta <- with(df,
  {
    if (method[2] == "both") {
      beta <- (Z + c(Z[-1], NA)) / (2 * S * I * 1)
    } else if (method[2] == "backward") {
      beta <- Z / (S * I * 1)
    } else if (method[2] == "forward") {
      beta <- c(Z[-1], NA) / (S * I * 1)
    } else {
      stop(
        "SI method: `method[2]` must be ",
        "\"forward\", \"backward\", or \"both\".",
        call. = FALSE
      )
    }
    beta
  }
)


## 7. Warn if `S` ever negative ----------------------------------------

if (any(df$S < 0, na.rm = TRUE)) {
  # May have underestimated `prep` or `nu`
  new_par_vals <- paste0(
    "\n* `prep` > ", sprintf("%.3f", prep),
    if (B_was_not_def) paste0("\n* `nu` > ", sprintf("%.3e", nu))
  )
  warning(
    "SI method: `S[i]` < 0 for some `i`. Retry with:",
    new_par_vals,
    call. = FALSE
  )
}


## 8. Warn if `I` is ever negative -------------------------------------

if (any(df$I < 0, na.rm = TRUE)) {
  # May have overestimated `prep` or `mu`, underestimated `tgen`
  new_par_vals <- paste0(
    "\n* `prep` < ", sprintf("%.3f", prep),
    if (mu_was_not_def) paste0("\n* `mu` < ", sprintf("%.3e", mu)),
    "\n* `tgen` > ", sprintf("%.1f", tgen)
  )
  warning(
    "SI method: `I[i]` < 0 for some `i`. Retry with:",
    new_par_vals,
    call. = FALSE
  )
}


attr(df, "par_list") <- par_list
attr(df, "method")   <- method
df
}
