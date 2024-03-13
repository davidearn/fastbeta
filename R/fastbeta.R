fastbeta <-
function (series, constants, ...)
{
	stopifnot(is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	          min(0, series, na.rm = TRUE) >= 0,
	          is.double(constants),
	          length(constants) == 5L,
	          is.finite(constants),
	          all(constants >= 0))
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	X <- .Call(R_fastbeta, series, constants)
	oldClass(X) <- oldClass(series)
	tsp(X) <- tsp(series)
	dimnames(X) <- list(NULL, c("S", "I", "R", "beta"))
	X
}
