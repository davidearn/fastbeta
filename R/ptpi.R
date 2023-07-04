ptpi <-
function (series, start, a = 0L, b = nrow(series) - 1L,
          tol = 1e-03, iter.max = 20L, complete = FALSE)
{
	stopifnot(exprs = {
		is.mts(series)
		is.double(series)
		ncol(series) == 3L
		min(0, series, na.rm = TRUE) >= 0
		is.numeric(start)
		length(start) == 1L
		start >= 0
		is.integer(a)
		length(a) == 1L
		a >= 0L
		is.integer(b)
		length(b) == 1L
		b < nrow(series)
		a < b
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
	storage.mode(start) <- "double"
	.Call(R_ptpi, series, start, a, b, tol, iter.max, complete)
}
