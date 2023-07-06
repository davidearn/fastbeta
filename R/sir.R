sir <-
function (n, beta, nu, mu, constants, stochastic = TRUE,
          prob = 1, delay = 1, useCompiled = TRUE, ...)
{
	stopifnot(exprs = {
		is.integer(n)
		length(n) == 1L
		n >= 1L
		n < .Machine$integer.max
		is.function(beta)
		!is.null(formals(beta))
		is.function(nu)
		!is.null(formals(nu))
		is.function(mu)
		!is.null(formals(mu))
		is.double(constants)
		length(constants) == 4L
		is.finite(constants)
		all(constants >= 0)
		is.double(prob)
		any(length(prob) == c(1L, n))
		min(prob) >= 0
		max(prob) <= 1
		is.double(delay)
		length(delay) >= 1L
		min(delay) >= 0
		sum(delay) >  0
	})

	if (stochastic) {

		tl.params <-
			(function (maxtau = 1, ...) list(maxtau = maxtau, ...))(...)

		init <- c(S = constants[[2L]],
		          I = constants[[3L]],
		          R = constants[[4L]],
		          B = 0,
		          Z = 0)
		tran <- list(c(S =  1, B = 1),        # birth
		             c(S = -1, I = 1, Z = 1), # infection
		             c(I = -1, R = 1),        # removal
		             c(S = -1),               # natural mortality
		             c(I = -1),               # ""
		             c(R = -1))               # ""
		if (useCompiled) {
			.Call(R_adsir_initialize, beta, nu, mu, constants[[1L]])
			on.exit(.Call(R_adsir_finalize), add = TRUE)
			ff <- function (x, theta, t) .Call(R_adsir_dot, t, x)
			Df <- function (x, theta, t) .Call(R_adsir_jac, t, x)
			X. <- ssa.adaptivetau(
				init.values  = init,
				transitions  = tran,
				rateFunc     = ff,
				params       = NULL,
				tf           = n,
				jacobianFunc = Df,
				tl.params    = tl.params)
		} else {
			Dm <- matrix(0, 4L, 6L)
			Di <- c(1L, 6L, 7L, 12L, 16L, 22L, 28L)
			ff <- function (x, theta, t) {
				xS <- x[[1L]]
				xI <- x[[2L]]
				xR <- x[[3L]]
				beta <- theta[[1L]](t)
				nu <- theta[[2L]](t)
				mu <- theta[[3L]](t)
				gamma <- theta[[4L]]
				c(nu,
				  beta * xS * xI,
				  if (xI > 1) gamma * xI else 0,
				  mu * xS,
				  if (xI > 1) mu * xI else 0,
				  mu * xR)
			}
			Df <- function (x, theta, t) {
				xS <- x[[1L]]
				xI <- x[[2L]]
				xR <- x[[3L]]
				beta <- theta[[1L]](t)
				nu <- theta[[2L]](t)
				mu <- theta[[3L]](t)
				gamma <- theta[[4L]]
				Dm[Di] <<- c(nu, beta * xI, beta * xS, gamma, mu, mu, mu)
				Dm
			}
			X. <- ssa.adaptivetau(
				init.values  = init,
				transitions  = tran,
				rateFunc     = ff,
				params       = list(beta, nu, mu, constants[[1L]]),
				tf           = n,
				jacobianFunc = Df,
				tl.params    = tl.params)
		}

		N <- dim(X.)[1L]
		i <- N - match(0L:n, as.integer(ceiling(X.[N:1L, 1L]))) + 1L
		if (anyNA(i)) {
			## MJ: maxtau constrains leaps but not steps => LOCF
			i[is.na(i)] <- 0L
			k <- c(which(i[-1L] != i[-length(i)]), length(i)) # run stops
			ik <- i[k]
			w <- which(ik == 0L)
			ik[w] <- ik[w - 1L]
			i <- rep.int(ik, k - c(0L, k[-length(k)]))
		}
		X <- X.[i, 2L:6L, drop = FALSE]

	} else {

		init <- c(S = constants[[2L]],
		          logI = log(constants[[3L]]),
		          R = constants[[4L]],
		          B = 0,
		          Z = 0)
		if (useCompiled) {
			.Call(R_desir_initialize, beta, nu, mu, constants[[1L]])
			on.exit(.Call(R_desir_finalize), add = TRUE)
			X. <- ode(
				y = init,
				times = 0:n,
				func = "R_desir_dot",
				parms = NULL,
				jacfunc = "R_desir_jac",
				jactype = "fullusr",
				hmax = 1,
				ynames = FALSE,
				dllname = "fastbeta",
				initfunc = NULL,
				initpar = NULL,
				...)
		} else {
			Dm <- matrix(0, 5L, 5L)
			Di <- c(1L, 2L, 5L, 6L, 8L, 10L, 13L)
			gg <- function (t, x, theta) {
				xS <- x[[1L]]
				xI <- exp(x[[2L]])
				xR <- x[[3L]]
				beta <- theta[[1L]](t)
				nu <- theta[[2L]](t)
				mu <- theta[[3L]](t)
				gamma <- theta[[4L]]
				beta.xS <- beta * xS
				beta.xS.xI <- beta.xS * xI
				list(c(nu - beta.xS.xI - mu * xS,
				       beta.xS - gamma - mu,
				       gamma * xI - mu * xR,
				       nu,
				       beta.xS.xI))
			}
			Dg <- function (t, x, theta) {
				xS <- x[[1L]]
				xI <- exp(x[[2L]])
				xR <- x[[3L]]
				beta <- theta[[1L]](t)
				nu <- theta[[2L]](t)
				mu <- theta[[3L]](t)
				gamma <- theta[[4L]]
				beta.xI <- beta * xI
				beta.xS.xI <- beta.xI * xS
				Dm[Di] <<- c(-beta.xI - mu,
				             beta,
				             beta.xI,
				             -beta.xS.xI,
				             gamma * xI,
				             beta.xS.xI,
				             -mu)
				Dm
			}
			X. <- ode(
				y = init,
				times = 0:n,
				func = gg,
				parms = list(beta, nu, mu, constants[[1L]]),
				jacfunc = Dg,
				jactype = "fullusr",
				hmax = 1,
				ynames = FALSE,
				dllname = NULL,
				initfunc = NULL,
				initpar = NULL,
				...)
		}

		X.[, 3L] <- exp(X.[, 3L])
		N <- dim(X.)[1L]
		X <-
			if (N < n + 1L)
				rbind(X.[, 2L:6L, drop = FALSE],
				      matrix(NA_real_, n + 1L - N, 5L))
			else X.[, 2L:6L, drop = FALSE]

	}
	head <- 1L:n
	tail <- 2L:(n + 1L)
	X[tail, 4L:5L] <- X[tail, 4L:5L] - X[head, 4L:5L]
	X[  1L, 4L:5L] <- NA_real_

	m.p <- missing(prob)
	m.d <- missing(delay)
	if (doObs <- !(m.p && m.d)) {
		X <- cbind(X, 0, deparse.level = 0L)
		Xt6 <- X[tail, 5L]
		if (stochastic) {
			Xt6 <- as.integer(Xt6)
			if (!m.p)
				Xt6 <- rbinom(n, Xt6, prob)
			if (!m.d)
				## FIXME? 'rmultinom' is more efficient, but not vectorized ...
				Xt6 <- tabulate(rep.int(1L:n, Xt6) +
				                sample(seq.int(0L, length.out = length(delay)),
				                       size = sum(Xt6),
				                       replace = TRUE,
				                       prob = delay),
				                n)
		} else {
			if (!m.p)
				Xt6 <- Xt6 * prob
			if (!m.d) {
				d <- length(delay) - 1L
				Xt6 <- filter(c(double(d), Xt6), delay / sum(delay),
				              sides = 1)[(d+1L):(d+n)]
			}
		}
		X[tail, 6L] <- Xt6
		X[  1L, 6L] <- NA_real_
	}
	ts(X, start = 0, names = c("S", "I", "R", "B", "Z", if (doObs) "Z.obs"))
}
