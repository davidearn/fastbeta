impute_zeros <- function(x, y) {

# Do nothing if there are no zeros in `y`
is_zero <- y == 0
if (!any(is_zero, na.rm = TRUE)) {
  return(y)
}

# Index zeros between first and last nonzero elements of `y`
index_zero <- which(is_zero)    
index_nonzero <- which(y != 0)
index_zero <- index_zero[
  index_zero > min(index_nonzero) &
  index_zero < max(index_nonzero)
]

# Linear interpolant of `(x, y)` with nonzero `y`
impute_y <- stats::approxfun(
  x      = x[index_nonzero],
  y      = y[index_nonzero],
  method = "linear",
  rule   = 1 # return `NA` outside the range of `x[index_nonzero]`
)

# `y` with zeros imputed by linear interpolation (where possible)
replace(y, index_zero, impute_y(x[index_zero]))
}

