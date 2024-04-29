seir <-
function (length.out = 1L, beta, nu, mu, constants, m = 0L, n = 1L,
          stochastic = TRUE, prob = 1, delay = 1, useCompiled = TRUE, ...)
{
	stopifnot(is.integer(m),
	          length(m) == 1L,
	          m >= 0L && m < 4096L,
	          is.integer(n),
	          length(n) == 1L,
	          n >= 1L && n < 4096L,
	          is.integer(length.out),
	          length(length.out) == 1L,
	          length.out >= 1L,
	          is.function(beta),
	          !is.null(formals(beta)),
	          is.function(nu),
	          !is.null(formals(nu)),
	          is.function(mu),
	          !is.null(formals(mu)),
	          is.double(constants),
	          length(constants) == m + n + 5L,
	          all(is.finite(constants)),
	          min(constants) >= 0,
	          is.double(prob),
	          any(length(prob) == c(1L, length.out - 1L)),
	          min(prob) >= 0,
	          max(prob) <= 1,
	          is.double(delay),
	          length(delay) >= 1L,
	          min(delay) >= 0,
	          sum(delay) > 0)

	p <- m + n + 2L

	init <- c(constants[seq_len(p)], 0, 0)
	names(init) <- nms <-
		c("S",
		  sprintf("E.%03x", seq_len(m)),
		  sprintf("I.%03x", seq_len(n)),
		  "R",
		  "Z", # incidence
		  "B") # births

	if (!useCompiled) {
		sigma <- constants[[p + 1L]] * m
		gamma <- constants[[p + 2L]] * n
		delta <- constants[[p + 3L]] * 1
		i.S <- 1L
		i.E <- seq.int(    2L, length.out = m)
		i.I <- seq.int(m + 2L, length.out = n)
		i.R <- p
	}

	if (stochastic) {

		tl.params <-
			(function (maxtau = 1, ...) list(maxtau = maxtau, ...))(...)

		tran <-
			c(list(`names<-`(c(-1, 1, 1), c("S", nms[2L], "Z"))),
			  list(`names<-`(c(1, 1), c("S", "B"))),
			  lapply(1L:p,
			         function(i) `names<-`(-1, nms[i])),
			  lapply(2L:(p - 1L),
			         function(i) `names<-`(c(-1, 1), nms[c(i, i + 1L)])),
			  list(`names<-`(c(1, -1), c("S", "R"))))
		## infection, birth, natural mortality, removal, loss of immunity

		if (useCompiled) {
			.Call(R_adseir_initialize, beta, nu, mu, constants, m, n)
			on.exit(.Call(R_adseir_finalize))
			ff <- function (x, theta, t) .Call(R_adseir_dot, t, x)
			Df <- function (x, theta, t) .Call(R_adseir_jac, t, x)
		}
		else {
			## D[i, j] = d(rate of transition j)/d(state i)
			##
			##     [S   ]  x  0  mu    0   0   0      0      0      0
			##     [E[i]]  0  0   0   mu   0   0  sigma      0      0
			##     [I[j]]  y  0   0    0  mu   0      0  gamma      0
			##     [R   ]  0  0   0    0   0  mu      0      0  delta
			##     [Z   ]  0  0   0    0   0   0      0      0      0
			##     [B   ]  0  0   0    0   0   0      0      0      0
			##
			## where x = beta * sum(I), y = beta * S
			D <- matrix(0, p + 2L, p + p + 1L)
			k.0 <- seq.int(from = m + 2L, by = 1L,
			               length.out = n)
			k.1 <- seq.int(from = (p + 2L) * 2L + 1L, by = p + 3L,
			               length.out = p)
			k.2 <- seq.int(from = (p + 2L) * (p + 2L) + 2L, by = p + 3L,
			               length.out = p - 1L)
			D[k.2] <- rep.int(c(sigma, gamma, delta), c(m, n, 1L))

			ff <-
			function (x, theta, t)
			{
				x.S <- x[i.S]
				x.E <- x[i.E]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				nu   <- nu  (t)
				mu   <- mu  (t)
				ok <- sum(x.E) + (s <- sum(x.I)) > 1
				c(        beta  * x.S * s           ,
				          nu                        ,
				          mu    * x.S               ,
				  if (ok) mu    * x.E else double(m),
				  if (ok) mu    * x.I else double(n),
				          mu    * x.R               ,
				          sigma * x.E               ,
				  if (ok) gamma * x.I else double(n),
				          delta * x.R               )
			}
			Df <-
			function (x, theta, t)
			{
				x.S <- x[i.S]
				x.I <- x[i.I]
				beta <- beta(t)
				mu   <- mu  (t)
				D[ 1L] <<- beta * sum(x.I)
				D[k.0] <<- beta * x.S
				D[k.1] <<- mu
				D
			}
		}

		X. <- ssa.adaptivetau(
			init.values  = init,
			transitions  = tran,
			rateFunc     = ff,
			params       = NULL,
			tf           = length.out - 1L,
			jacobianFunc = Df,
			tl.params    = tl.params)

		D. <- dim(X.)
		i <- D.[1L] - match(seq.int(0L, length.out = length.out),
		                    as.integer(ceiling(X.[D.[1L]:1L, 1L]))) + 1L
		if (anyNA(i)) {
			## tl.params[["maxtau"]] constrains leaps but not steps => LOCF
			i[is.na(i)] <- 0L
			k <- c(which(i[-1L] != i[-length(i)]), length(i)) # run stops
			ik <- i[k]
			w <- which(ik == 0L)
			ik[w] <- ik[w - 1L]
			i <- rep.int(ik, k - c(0L, k[-length(k)]))
		}
		X <- X.[i, -1L, drop = FALSE] # discarding time

	}
	else {

		## E[i], I[j], R: logarithmic scale
		stopifnot(min(constants[2L:(m + n + 2L)]) > 0)
		init[2L:p] <- log(init[2L:p])

		if (useCompiled) {
			.Call(R_deseir_initialize, beta, nu, mu, constants, m, n)
			on.exit(.Call(R_deseir_finalize))
			gg <- "R_deseir_dot"
			Dg <- "R_deseir_jac"
		}
		else {
			## D[i, j] = d(rate of change in state i)/d(state j)
			##
			## nonzero pattern in m = 4, n = 4 case :
			##
			##       S E E E E I I I I R Z B
			##     S | . . . . | | | | | . .
			##     E | | . . . | | | | . . .
			##     E . | | . . . . . . . . .
			##     E . . | | . . . . . . . .
			##     E . . . | | . . . . . . .
			##     I . . . . | | . . . . . .
			##     I . . . . . | | . . . . .
			##     I . . . . . . | | . . . .
			##     I . . . . . . . | | . . .
			##     R . . . . . . . . | | . .
			##     Z | . . . . | | | | . . .
			##     B . . . . . . . . . . . .
			##
			## nonzero pattern in m = 0, n = 4 case :
			##
			##       S I I I I R Z B
			##     S | | | | | | . .
			##     I | . | | | . . .
			##     I . | | . . . . .
			##     I . . | | . . . .
			##     I . . . | | . . .
			##     R . . . . | | . .
			##     Z | | | | | . . .
			##     B . . . . . . . .
			##
			## where log(.) is suppressed only for pretty printing
			D <- matrix(0, p + 2L, p + 2L)
			k.0 <- seq.int(from = (p + 2L) * 2L + 3L, by = p + 3L,
			               length.out = p - 2L)
			k.1 <- k.0 - p - 2L

			i.1 <- 2L:(p - 1L)
			i.0 <- 3L:p
			a.1 <- rep.int(c(sigma, gamma), c(m, n))
			a.0 <- c(a.0[-1L], delta)

			gg <-
			function (t, x, theta)
			{
				x.S <- x[i.S]
				x.E <- x[i.E]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				nu   <- nu  (t)
				mu   <- mu  (t)
				s.0 <- sum(exp(x.I        ))
				s.1 <- sum(exp(x.I - x[2L]))
				list(c(nu - beta * x.S * s.0 + delta * exp(x.R) - mu * x.S,
				       beta * x.S * s.1 - a.1[1L] - mu,
				       a.1 * exp(x[i.1] - x[i.0]) - a.0 - mu,
				       beta * x.S * s.0,
				       nu))
			}
			Dg <-
			function (t, x, theta)
			{
				x.S <- x[i.S]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				mu   <- mu  (t)
				s.0 <- sum(u.0 <- exp(x.I        ))
				s.1 <- sum(u.1 <- exp(x.I - x[2L]))
				D[i.S, i.S] <<- -(D[p + 1L, i.S] <<- beta * s.0) - mu
				D[i.S, i.I] <<- -(D[p + 1L, i.I] <<- beta * x.S * u.0)
				D[i.S, i.R] <<- delta * exp(x.R)
				D[ 2L, i.S] <<- beta * s.1
				D[ 2L, i.I] <<- beta * x.S * u.1
				D[ 2L,  2L] <<- if (m) -beta * s.1 else 0
				D[k.0] <<- -(D[k.1] <<- a.1 * exp(x[i.1] - x[i.0]))
				D
			}
		}

		X. <- lsoda(
			y        = init,
			times    = seq.int(0, length.out = length.out),
			func     = gg,
			parms    = NULL,
			jacfunc  = Dg,
			jactype  = "fullusr",
			hmax     = 1,
			ynames   = FALSE,
			dllname  = if (useCompiled) "fastbeta",
			initfunc = NULL,
			initpar  = NULL,
			...)

		X. <- X.[, -1L, drop = FALSE]
		X.[, 2L:p] <- exp(X.[, 2L:p])
		D. <- dim(X.)
		X <-
			if (D.[1L] < length.out)
				## not seen, but help("lsoda") says that it is possible
				rbind(X., matrix(NaN, length.out - D.[1L], D.[2L]),
				      deparse.level = 0L)
			else X.

	}

	if (length.out > 1L) {
	head <- 1L:(length.out - 1L)
	tail <- 2L:length.out
	X[tail, p + 1L:2L] <- X[tail, p + 1L:2L] - X[head, p + 1L:2L]
	}
	X[  1L, p + 1L:2L] <- NA_real_

	m.p <- missing(prob)
	m.d <- missing(delay)
	if (doObs <- !(m.p && m.d)) {
		X <- cbind(X, NA_real_, deparse.level = 0L)
		if (length.out > 1L) {
		Z <- X[tail, p + 1L]
		if (stochastic) {
			Z <- as.integer(Z)
			if (!m.p)
				Z <- rbinom(length.out - 1L, Z, prob)
			if (!m.d)
				## FIXME? 'rmultinom' is more efficient, but not vectorized ...
				Z <- tabulate(rep.int(seq_len(length.out - 1L), Z) +
				              sample(seq.int(0L, length.out = length(delay)),
				                     size = sum(Z),
				                     replace = TRUE,
				                     prob = delay),
				              length.out - 1L)
		}
		else {
			if (!m.p)
				Z <- Z * prob
			if (!m.d) {
				d <- length(delay) - 1L
				Z <- filter(c(double(d), Z), delay / sum(delay),
				            sides = 1)[(d + 1L):(d + length.out - 1L)]
			}
		}
		X[tail, p + 3L] <- Z
		}
	}

	nms <- rep.int(c("S", "E", "I", "R", "Z", "B", "Z.obs"),
	               c(1L, m, n, 1L, 1L, 1L, if (doObs) 1L else 0L))
	ts(X, start = 0, names = nms)
}
