ptpi <- function(df = data.frame(), par_list = list(),
                 a, b,
                 initial_S0_est,
                 iter = 0L) {

## 1. Set-up -----------------------------------------------------------

# Assume constant vital rates if vital data were not supplied
if (is.null(df$B)) {
  df$B  <- with(par_list, nu * hatN0 * 1)
}
if (is.null(df$mu)) {
  df$mu <- with(par_list, mu)
}

# Preallocate memory for all susceptible time series,
# and initialize the first
S_mat <- matrix(NA, nrow = nrow(df), ncol = iter + 1)
S_mat[1, 1] <- initial_S0_est
S_mat[a, 1] <- initial_S0_est


## 2. Peak-to-peak iteration -------------------------------------------

for (j in seq_len(iter + 1)) {

  ## 2.(a) Update `SA` estimate

  # Reconstruct from index `a` to end
  for (i in (a+1):nrow(S_mat)) {
    S_mat[i, j] <- with(df[c("Z", "B", "mu")],
      {
        ((1 - 0.5 * mu[i-1] * 1) * S_mat[i-1,j] + B[i] - Z[i]) /
          (1 + 0.5 * mu[i] * 1)
      }
    )
  }
  if (j == iter + 1) {
    break
  }
  S_mat[a, j+1] <- S_mat[b, j]

  ## 2.(b) Update `S0` estimate

  # Reconstruct from index `a` to start (backwards in time)
  for (i in (a-1):1) {
    S_mat[i, j+1] <- with(df[c("Z", "B", "mu")],
      {
        ((1 + 0.5 * mu[i+1] * 1) * S_mat[i+1, j+1] - B[i+1] + Z[i+1]) /
          (1 - 0.5 * mu[i] * 1)
      }
    )
  }

}


out <- list(
  S_mat     = S_mat,
  S0       = S_mat[1, ],
  S0_final = S_mat[1, ncol(S_mat)],
  SA       = S_mat[a, ],
  SA_final = S_mat[a, ncol(S_mat)]
)
attr(out, "arg_list") <- as.list(environment())[names(formals(ptpi))]
out
}
