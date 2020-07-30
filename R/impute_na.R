impute_na <- function(x, y) {
  # Do nothing if there are no `NA` in `y`
  is_na <- is.na(y)
  if (!any(is_na)) {
    return(y)
  }

  # `y` with `NA` imputed by linear interpolation (where possible)
  stats::approx(
    x      = x[!is_na],
    y      = y[!is_na],
    xout   = x,
    method = "linear",
    rule   = 1 # return `NA` outside the range of `x[!is_na]`
  )
}
