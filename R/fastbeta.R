fastbeta <-
function (Z, B, mu, gamma, S0, I0)
{
	stopifnot(exprs = {
		is.numeric(Z)
		min(0, Z, na.rm = TRUE) >= 0
		is.numeric(B)
		length(B) == length(Z)
		min(0, B, na.rm = TRUE) >= 0
		is.double(mu)
		length(mu) == length(Z)
		min(0, mu, na.rm = TRUE) >= 0
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
	storage.mode(Z) <- storage.mode(B) <-
		storage.mode(S0) <- storage.mode(I0) <- "double"
	ts(.Call(R_fastbeta, Z, B, mu, gamma, S0, I0),
       start = 0, names = c("S", "I", "beta"))
}
