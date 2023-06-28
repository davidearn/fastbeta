fastbeta <-
function (series, constants)
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
	X <- .Call(R_fastbeta, series, constants)
	oldClass(X) <- oldClass(series)
	tsp(X) <- tsp(series)
	dimnames(X) <- list(NULL, c("beta", "S", "I"))
	X
}
