#' Compute RRMSE in a vector
#'
#' `compute_rrmse()` computes the relative root mean square error
#' (RRMSE) in a vector.
#'
#' @details
#' Let `cc` ("complete cases") be the vector of indices `i` such that
#' neither `x[i]` nor `xhat[i]` is missing (`NA`, `NaN`, `Inf`). If
#' `length(cc) = 0` or if `na_rm = FALSE` and `length(cc) < length(x)`,
#' then `compute_rrmse()` returns `NA`. Otherwise, it returns
#' `sqrt(sum((xhat[cc] - x[cc])^2) / length(cc)) / mean(x[cc])`.
#' This is the RRMSE in `xhat[cc]` unless `mean(x[cc]) = 0`,
#' in which case RRMSE is not defined and the expression is
#' (correctly) evaluated as `Inf` or `NaN`.
#' 
#' @param x Numeric vector.
#' @param xhat Numeric vector with length equal to `length(x)`.
#'   An estimate of `x`.
#' @param na_rm Logical. If `TRUE`, then missing values
#'   (`NA`, `NaN`, `Inf`) are ignored.
#'
#' @return
#' A numeric scalar indicating the RRMSE in `xhat`, with exceptions
#' (see Details).
#'
#' @examples
#' x <- rep(10, 100)
#' xhat <- x + rnorm(x, 0, 1)
#' compute_rrmse(x, xhat) # RRMSE in `xhat`
#' 
#' @md
#' @export
compute_rrmse <- function(x, xhat, na_rm = FALSE) {

cc <- !is.na(x) & !is.na(xhat) & !is.infinite(x) & !is.infinite(xhat)
if (sum(cc) == 0 || (sum(cc) < nrow(df) && !na_rm)) {
  return(NA)
}
mse <- sum((xhat[cc] - x[cc])^2) / sum(cc)
rrmse <- sqrt(mse) / mean(x[cc])
rrmse
}
