ptpi <-
function (series, constants, a = 0L, b = nrow(series) - 1L,
          tol = 1e-03, iter.max = 32L, complete = FALSE, ...)
{
	stopifnot(exprs = {
		is.mts(series)
		is.double(series)
		ncol(series) == 3L
		min(0, series, na.rm = TRUE) >= 0
		is.double(constants)
		length(constants) == 5L
		is.finite(constants)
		all(constants >= 0)
		is.numeric(a)
		length(a) == 1L
		a >= tsp(series)[1L]
		is.numeric(b)
		length(b) == 1L
		b <= tsp(series)[2L]
		b - a >= 1 / tsp(series)[3L]
		is.double(tol)
		length(tol) == 1L
		!is.na(tol)
		is.integer(iter.max)
		length(iter.max) == 1L
		iter.max >= 1L
		is.logical(complete)
		length(complete) == 1L
		!is.na(complete)
	})
	tsp <- tsp(series)
	a <- as.integer(round((a - tsp[1L]) * tsp[3L]))
	b <- as.integer(round((b - tsp[1L]) * tsp[3L]))
	if (a > b)
		stop("'a' is greater than 'b' after rounding; should never happen ...")
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	r <- .Call(R_ptpi, series, constants, a, b, tol, iter.max, complete)
	if (complete) {
		oldClass(r[["X"]]) <- oldClass(series)
		tsp(r[["X"]]) <- c(tsp[1L] + c(a, b) / tsp[3L], tsp[3L])
	}
	r
}
