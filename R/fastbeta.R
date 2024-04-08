fastbeta <-
function (series, constants, m = 0L, n = 1L, ...)
{
	stopifnot(is.integer(m),
	          length(m) == 1L,
	          m >= 0L,
	          is.integer(n),
	          length(n) == 1L,
	          n >= 1L,
	          is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	          min(0, series, na.rm = TRUE) >= 0,
	          is.double(constants),
	          length(constants) == m + n + 5L,
	          all(is.finite(constants)),
	          min(constants) >= 0)
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	X <- .Call(R_fastbeta, series, constants, m, n)
	oldClass(X) <- oldClass(series)
	tsp(X) <- tsp(series)
	dimnames(X) <-
		list(NULL, rep.int(c("S", "E", "I", "R", "beta"), c(1L, m, n, 1L, 1L)))
	X
}
