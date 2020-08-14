#' Impute missing values in a numeric vector
#'
#' @description
#' Imputes missing values in a numeric vector
#' by linear interpolation of observed values,
#' assuming equal spacing (e.g., that the
#' elements of the vector are observations
#' of a quantity at equally spaced times).
#'
#' @param x Numeric vector.
#' @param zero_as_na Logical scalar. If `TRUE`,
#'   then zeros in `x` are treated like missing values.
#'
#' @return
#' `x` with missing values replaced where possible.
#'
#' @details
#' Missing values include `NA`, `NaN`, `Inf`, and `-Inf`.
#' Those not preceded by or not followed by an observed
#' value cannot be imputed by linear interpolation and
#' are replaced with the common value `NA`.
#'
#' If `zero_as_na` is `TRUE`, then zeros preceded by and
#' followed by an observed value (zero or nonzero) are
#' treated like `NA`. Zeros not preceded by or not followed
#' by an observed value (zero or nonzero) are left untouched.
#'
#' @examples
#' x <- c(NA, 0, NA, 0, NA, 2, NA, 4, 0, 6, NaN, Inf, -Inf)
#' impute_na(x)
#' impute_na(x, zero_as_na = TRUE)
#'
#' @export
impute_na <- function(x, zero_as_na = FALSE) {
  # Do nothing if there is nothing to do
  if (!is.numeric(x)) {
    return(x)
  }
  if (all(is.finite(x))) {
    if (!zero_as_na || (zero_as_na && !any(x == 0))) {
      return(x)
    }
  }

  x[is.nan(x) | is.infinite(x)] <- NA

  # Do nothing else if nothing else can be done
  r <- range(which(!is.na(x)))
  if (r[2] - r[1] < 2) {
    return(x)
  }

  y <- x[r[1]:r[2]]
  if (zero_as_na) {
    y_end <- y[c(1, length(y))]
    y[y == 0] <- NA
    y[c(1, length(y))] <- y_end
  }
  approx_out <- stats::approx(
    x      = which(!is.na(y)),
    y      = y[!is.na(y)],
    xout   = seq_along(y),
    method = "linear"
  )
  y <- approx_out$y
  x[r[1]:r[2]] <- y
  x
}
