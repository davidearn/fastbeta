bootstrap <-
function (r, series, constants, ...)
{
	stopifnot(exprs = {
		is.mts(series)
		is.double(series)
		ncol(series) == 3L
		min(0, series, na.rm = TRUE) >= 0
		is.double(constants)
		length(constants) == 3L
		is.finite(constants)
		all(constants >= 0)
	})

	nu. <-
	    if(TRUE)
	        suppressWarnings(
	            c(1, -1) * cumsum(c(1, rep_len(c(-2, 2), n)) * series[, 2L]))
	    else {
	        L <- suppressWarnings(
	            matrix(c(0.5, 0.5, double(n)), n + 1L, n + 1L))
	        L[1L] <- 1
	        forwardsolve(L, series[, 2L])
	    }

	n <- nrow(series) - 1L
	s <- as.double(0:n)
	beta <- approxfun(s, series[, 1L],
	                  method = "linear", rule = 2:1, ties = "ordered")
	nu   <- approxfun(s, nu.,
	                  method = "linear", rule = 2:1, ties = "ordered")
	mu   <- approxfun(s, series[, 3L],
	                  method = "linear", rule = 2:1, ties = "ordered")

	constants. <- c(constants, 0) # R0 is irrelevant
	R <- replicate(r, {
		X <- sir(n, beta, nu, mu, constants., ...)
		series[, 1L:2L] <<- X[, c(4L, ncol(X))]
		fastbeta(series, constants)[, 1L]
	})
	oldClass(R) <- oldClass(series)
	tsp(R) <- tsp(series)
	R
}
