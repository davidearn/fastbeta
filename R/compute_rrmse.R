compute_rrmse <- function(x, y, na_rm = TRUE) {

df <- data.frame(x, y)
cc <- apply(df, 1, function(x) !any(is.na(x) | is.infinite(x)))

if (sum(cc) == 0 || (sum(cc) < nrow(df) && !na_rm)) {
  return(NA)
}

df <- df[cc, ]
with(df,
  {
    mse <- sum((x - y)^2) / length(x)
    rrmse <- sqrt(mse) / mean(x)
    rrmse
  }
)
}
