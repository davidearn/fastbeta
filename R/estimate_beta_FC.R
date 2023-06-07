estimate_beta_FC <- function(df       = data.frame(),
                             par_list = list()) {

## 1. Set-up -----------------------------------------------------------

# Load necessary elements of `par_list` into the execution environment
list2env(
  par_list[c("prep", "trep", "S0", "hatN0", "nu", "tgen")],
  envir = environment()
)

# Logical: Is `df` missing a column `B`?
B_was_not_def <- is.null(df$B)

# Preallocate memory for output. Assume a constant birth rate,
# if necessary.
df <- data.frame(
  t     = df$t,
  C     = df$C,
  Z     = NA,
  Z_agg = NA,
  B     = if (B_was_not_def) nu * hatN0 * 1 else df$B,
  B_agg = NA,
  S     = NA,
  I     = NA,
  beta  = NA
)


## 2. Impute missing values in reported incidence, births --------------

df[c("C", "B")] <- lapply(df[c("C", "B")],
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
##       in `df$C`. In the worst case, `df$C` is now looks like
##       `c(NA,...,NA,0,...,0,+,...,+,0,...,0,NA,...,NA)`,
##       where `+` denotes a positive number.


## 4. Estimate incidence -----------------------------------------------

trepr <- round(trep)
df$Z <- c(
  (1 / prep) * df$C[(trepr+1):nrow(df)],
  rep(NA, trepr)
)


## 5. Aggregate incidence, births --------------------------------------

# Aggregates are reported after each mean generation interval,
# starting at the first time point. Aggregates at the first
# time point cannot be computed, because there are no prior
# observations to aggregate.
tgenr <- round(tgen)
ind_with_agg <- seq(1, nrow(df), by = tgenr)
for (i in ind_with_agg[-1]) {
  ind_into_agg <- seq(i - tgenr + 1, i, by = 1)
  df$B_agg[i] <- sum(df$B[ind_into_agg])
  df$Z_agg[i] <- sum(df$Z[ind_into_agg])
}


## 6. Estimate susceptibles, infecteds, transmission rate --------------

df[c("S", "I", "beta")] <- with(df,
  {
    S[1] <- S0
    for (i in ind_with_agg[-1]) {
      S[i] <- S[i-tgenr] + B_agg[i] - Z_agg[i]
    }
    I <- Z_agg
    beta[ind_with_agg] <- c(Z_agg[ind_with_agg[-1]], NA) /
      (S[ind_with_agg] * Z_agg[ind_with_agg] * tgenr)
    list(S, I, beta)
  }
)


## NOTE: If too many `NA` were retained at the start of `df$C` or
##       `df$B`, then `df$S` and `df$beta` will be filled with `NA`.
##       If zeros were retained in `df$C`, then `NaN` or `Inf` will
##       appear in `df$beta` whenever a divide-by-zero error occurs.


## 7. Warn if `S` is ever negative -------------------------------------

if (any(df$S < 0, na.rm = TRUE)) {
  # May have underestimated `prep` or `nu`
  new_par_vals <- paste0(
    "\n* `prep` > ", sprintf("%.3f", prep),
    if (B_was_not_def) paste0("\n* `nu` > ", sprintf("%.3e", nu))
  )
  warning(
    "FC method: `S[i]` < 0 for some `i`. Retry with:",
    new_par_vals,
    call. = FALSE
  )
}


## 8. Warn if `beta` is ever `NaN` or `Inf` ----------------------------

if (any(is.nan(df$beta) | is.infinite(df$beta))) {
  warning(
    "FC method: `beta[i]` is `NaN` or `Inf` for some `i`. \n",
    "Is the first or last observation in `df$C` a zero?",
    call. = FALSE
  )
}


attr(df, "par_list") <- par_list
df
}
