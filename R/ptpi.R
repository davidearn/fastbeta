ptpi <-
function (series, sigma = gamma, gamma = 1, delta = 0, init,
          m = length(init) - n - 2L, n = 1L, a = 0L, b = nrow(series) - 1L,
          tol = 1e-03, iter.max = 32L,
          complete = FALSE, backcalc = FALSE, ...)
{
	stopifnot(is.integer(m) && length(m) == 1L && m >= 0L,
	          is.integer(n) && length(n) == 1L && n >= 1L,
	          is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	          min(0, series, na.rm = TRUE) >= 0,
	          is.double(sigma) && length(sigma) == 1L && sigma >= 0,
	          is.double(gamma) && length(gamma) == 1L && gamma >= 0,
	          is.double(delta) && length(delta) == 1L && sigma >= 0,
	          is.double(init),
	          length(init) == m + n + 2L,
	          all(is.finite(init)),
	          min(init) >= 0,
	          is.numeric(a),
	          length(a) == 1L,
	          a >= tsp(series)[1L],
	          is.numeric(b),
	          length(b) == 1L,
	          b <= tsp(series)[2L],
	          b - a >= 1 / tsp(series)[3L],
	          is.double(tol),
	          length(tol) == 1L,
	          !is.na(tol),
	          is.integer(iter.max),
	          length(iter.max) == 1L,
	          iter.max >= 1L,
	          is.logical(complete),
	          length(complete) == 1L,
	          !is.na(complete),
	          is.logical(backcalc),
	          length(backcalc) == 1L,
	          !is.na(backcalc))
	tsp <- tsp(series)
	a <- as.integer(round((a - tsp[1L]) * tsp[3L]))
	b <- as.integer(round((b - tsp[1L]) * tsp[3L]))
	if (a > b)
		stop(gettextf("'%s' is greater than '%s' after rounding; should never happen ...",
		              "a", "b"),
		     domain = NA)
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	r <- .Call(R_ptpi, series, sigma, gamma, delta, init,
	           m, n, a, b, tol, iter.max, complete, backcalc)
	names(r[["value"]]) <- rep.int(c("S", "E", "I", "R"), c(1L, m, n, 1L))
	if (complete) {
		oldClass(r[["X"]]) <- oldClass(series)
		tsp(r[["X"]]) <- c(tsp[1L] + c(a, b) / tsp[3L], tsp[3L])
		dimnames(r[["X"]]) <- list(NULL, names(r[["value"]]), NULL)
	}
	r
}
