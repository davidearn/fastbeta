ptpi <-
function (data = ts(cbind(Z, B, mu), start = 0),
          Z, B, mu, start, a = 0L, b = nrow(data) - 1L,
          tol = 1e-06, iter.max = 20L, complete = FALSE)
{
	stopifnot(exprs = {
		is.mts(data)
		is.double(data)
		ncol(data) == 3L
		min(0, data, na.rm = TRUE) >= 0
		is.numeric(start)
		length(start) == 1L
		start >= 0
		is.integer(a)
		length(a) == 1L
		a >= 0L
		is.integer(b)
		length(b) == 1L
		b < nrow(data)
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
	.Call(R_ptpi, data, start, a, b, tol, iter.max, complete)
}
