fastbeta.affine <-
function (series, constants, m, n)
{
	r <- nrow(series)
	p <- m + n + 2L

	Z    <- series[-1L, 1L]
	B    <- series[-1L, 2L]
	mu.0 <- series[-1L, 3L]
	mu.1 <- series[- r, 3L]

	sigma <- constants[p + 1L]
	gamma <- constants[p + 2L]
	delta <- constants[p + 3L]

	hsigma <- 0.5 * m * sigma
	hgamma <- 0.5 * n * gamma
	hdelta <- 0.5     * delta
	hmu.0  <- 0.5     * mu.0
	hmu.1  <- 0.5     * mu.1

	hL <- rep.int(c(   hsigma, hgamma, hdelta, 0), c(    m, n, 1L, 1L))
	hR <- rep.int(c(1, hsigma, hgamma, hdelta   ), c(1L, m, n, 1L    ))

	## A := [A0, A1]
	## y := [ 1;  x]
	## x := A * y = A0 + A1 * x

	make.A <-
	function (j)
	{
		A <- matrix(0, p, p + 1L)

		rL <- (1 - hL - hmu.0[j]) / (1 + hL + hmu.1[j])
		rR <- hR / (1 + hL + hmu.1[j])

		## set band(A,  0,  0)
		A[seq.int(from = p + 1L, by = p + 1L, length.out = p     )] <- rL
		## set band(A, -1, -1)
		A[seq.int(from = p + 2L, by = p + 1L, length.out = p - 1L)] <- rR[2L:p]
		## set b[1]
		A[1L, 1L] <- rR[1L] * Z[j]
		## set b[p]
		A[ p, 1L] <- (B[j] - Z[j]) / (1 + hmu.1[j])

		tmp <- A[1L, ]
		for (i in 2L:p)
			A[i, ] <- tmp <- A[i, ] + rR[i] * tmp

		A
	}

	l.A <- vector("list", r - 1L)
	l.x <- vector("list", r)
	y <- rep.int(1, p + 1L)

	l.x[[1L]] <- y[-1L] <- constants[c(2L:p, 1L)]
	for (j in seq_len(r - 1L))
		l.x[[j + 1L]] <- y[-1L] <- as.double((l.A[[j]] <- make.A(j)) %*% y)

	beta <- 0.5 * (series[-r, 1L] + Z) /
		vapply(l.x[-1L], function(x) x[p] * sum(x[(p - 1L - 1L):(p - 1L - n)]), 0)

	list(beta, l.x, l.A, l.x)
}
