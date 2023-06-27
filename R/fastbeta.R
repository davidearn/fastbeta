fastbeta <-
function (data = ts(cbind(Z, B, mu), start = 0),
          Z, B, mu, gamma, S0, I0)
{
	stopifnot(exprs = {
		is.mts(data)
		is.double(data)
		ncol(data) == 3L
		min(0, data, na.rm = TRUE) >= 0
		is.double(gamma)
		length(gamma) == 1L
		gamma >= 0
		is.numeric(S0)
		length(S0) == 1L
		S0 >= 0
		is.numeric(I0)
		length(I0) == 1L
		I0 >= 0
	})
	storage.mode(S0) <- storage.mode(I0) <- "double"
	X <- .Call(R_fastbeta, data, gamma, S0, I0)
	oldClass(X) <- oldClass(data)
	tsp(X) <- tsp(data)
	dimnames(X) <- list(NULL, c("S", "I", "beta"))
	X
}
