#' Compute RRMSE in a time series estimate
#'
#' `compute_rrmse()` computes the relative root mean square error
#' (RRMSE) in a time series estimate of a known time-varying quantity.
#'
#' @param act Numeric vector. A time series.
#' @param est Numeric vector with length equal to `length(act)`.
#'   A time series estimate of `act`.
#' @param na_rm Logical. If `TRUE`, then missing values
#'   (`NA`, `NaN`, `Inf`) are ignored.
#'
#' @return
#' Let `cc` ("complete cases") be the vector of indices `i` such that
#' neither `act[i]` nor `est[i]` is missing. If `length(cc) = 0` or if
#' `length(cc) < length(act)` and `na_rm = FALSE`, then `NA`. Otherwise,
#' the RRMSE in `est` given by
#' `sqrt(sum(act[cc] - est[cc])^2 / length(cc)) / mean(act[cc])`.
#'
#' @md
#' @export
compute_rrmse <- function(act, est, na_rm = TRUE) {

df <- data.frame(act, est)
is_missing <- apply(df, 1, function(x) any(is.na(x) | is.infinite(x)))
cc <- which(!is_missing) # complete cases

if (length(cc) == 0 || (length(cc) < nrow(df) && !na_rm)) {
  return(NA)
}

df <- df[cc, ]
with(df,
  {
    mse <- sum((act - est)^2) / nrow(df)
    rrmse <- sqrt(mse) / mean(act)
    rrmse
  }
)
}
