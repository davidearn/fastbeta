#' @title Compute RRMSE in a vector
#'
#' @description
#' Computes the relative root mean square error (RRMSE) in a vector.
#' 
#' @param x Numeric vector.
#' @param xhat Numeric vector with length equal to `length(x)`.
#'   An estimate of `x`.
#' @param na_rm Logical. If `TRUE`, then missing values are ignored.
#'
#' @return
#' A numeric scalar indicating the RRMSE in `xhat`, with exceptions
#' (see Details).
#'
#' @details
#' Let `cc = is.finite(x) & is.finite(xhat)` ("complete cases"). `cc` is the
#' vector of indices `i` such that both `x[i]` and `xhat[i]` are finite
#' (not missing and not infinite). If `length(cc) = 0` or if `na_rm = FALSE`
#' and `length(cc) < length(x)`, then `compute_rrmse()` returns `NA`. Otherwise,
#' it returns `sqrt(sum((xhat[cc] - x[cc])^2) / length(cc)) / mean(x[cc])`.
#' This is the RRMSE in `xhat[cc]` unless `mean(x[cc]) = 0`, in which case
#' RRMSE is not defined and the expression is evaluated accordingly as `Inf`
#' or `NaN`.
#'
#' @examples
#' x <- rep(10, 100)
#' xhat <- x + rnorm(x, 0, 1)
#' compute_rrmse(x, xhat) # RRMSE in `xhat`
#' 
#' @export
compute_rrmse <- function(x, xhat, na_rm = FALSE) {
  cc <- is.finite(x) & is.finite(xhat)
  if (sum(cc) == 0 || (sum(cc) < length(x) && !na_rm)) {
    return(NA)
  }
  mse <- sum((xhat[cc] - x[cc])^2) / sum(cc)
  sqrt(mse) / mean(x[cc])
}
