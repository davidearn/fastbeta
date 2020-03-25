#' Compute RRMSE in a time series estimate
#'
#' `compute_rrmse()` computes the relative root mean square error
#' (RRMSE) in a time series estimate of a known time-varying quantity.
#'
#' @param x Numeric vector. A time series.
#' @param y Numeric vector with length equal to `length(x)`.
#'   A time series estimate of `x`.
#' @param na_rm Logical. If `TRUE`, then missing values
#'   (`NA`, `NaN`, `Inf`) are ignored.
#'
#' @return
#' Let `cc` ("complete cases") be the vector of indices `i` such that
#' neither `x[i]` nor `y[i]` is missing. If `length(cc) = 0` or if
#' `length(cc) < length(x)` and `na_rm = FALSE`, then `NA`. Otherwise,
#' the RRMSE in `y` given by
#' `sqrt(sum((x[cc] - y[cc])^2) / length(cc)) / mean(x[cc])`.
#'
#' @examples
#' x <- rep(c(0, 1, 0, -1), 100)
#' y <- x + rnorm(x, 0, 0.01)
#' compute_rrmse(x, y) # RRMSE in `y`
#' 
#' @md
#' @export
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
