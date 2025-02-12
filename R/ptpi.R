ptpi <-
function (series, sigma = 1, gamma = 1, delta = 0,
          m = 1L, n = 1L, init,
          start = tsp(series)[1L], end = tsp(series)[2L],
          tol = 1e-03, iter.max = 32L,
          backcalc = FALSE, complete = FALSE, ...)
{
	stopifnot(is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	          min(0, series, na.rm = TRUE) >= 0,
	          is.double(sigma) && length(sigma) == 1L && sigma >= 0,
	          is.double(gamma) && length(gamma) == 1L && gamma >= 0,
	          is.double(delta) && length(delta) == 1L && sigma >= 0,
	          is.integer(m) && length(m) == 1L && m >= 0L,
	          is.integer(n) && length(n) == 1L && n >= 1L,
	          is.double(init),
	          length(init) == m + n + 2L,
	          all(is.finite(init)),
	          min(init) >= 0,
	          is.numeric(start),
	          length(start) == 1L,
	          start >= tsp(series)[1L],
	          is.numeric(end),
	          length(end) == 1L,
	          end < tsp(series)[2L] + 1 / tsp(series)[3L],
	          start <= end,
	          is.double(tol),
	          length(tol) == 1L,
	          !is.na(tol),
	          is.integer(iter.max),
	          length(iter.max) == 1L,
	          iter.max >= 1L,
	          is.logical(backcalc),
	          length(backcalc) == 1L,
	          !is.na(backcalc),
	          is.logical(complete),
	          length(complete) == 1L,
	          !is.na(complete))
	tsp <- tsp(series)
	a <- as.integer(trunc((start - tsp[1L]) * tsp[3L]))
	b <- as.integer(trunc((  end - tsp[1L]) * tsp[3L])) + 1L
	if (a >= b)
		stop(gettextf("'%s' is greater than '%s' after truncation; should never happen ...",
		              "start", "end"),
		     domain = NA)
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	ans <- .Call(R_ptpi, series, sigma, gamma, delta, m, n, init,
	             a, b, tol, iter.max, backcalc, complete)
	names(ans[["value"]]) <- rep(c("S", "E", "I", "R"), c(1L, m, n, 1L))
	if (complete) {
		oldClass(ans[["x"]]) <- c("mts", "ts", "array")
		tsp(ans[["x"]]) <- c(tsp[1L] + c(a, b - 1L) / tsp[3L], tsp[3L])
		dimnames(ans[["x"]]) <- list(NULL, names(ans[["value"]]), NULL)
	}
	ans
}
