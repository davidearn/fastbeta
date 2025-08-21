fastbeta <-
function (series, sigma = 1, gamma = 1, delta = 0,
          m = 1L, n = 1L, init, ...)
{
	stopifnot(is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	       ## min(0, series         , na.rm = TRUE) >= 0, # original
	          min(0, series[, 1L:2L], na.rm = TRUE) >= 0, # experiment
	          is.double(sigma) && length(sigma) == 1L && sigma >= 0,
	          is.double(gamma) && length(gamma) == 1L && gamma >= 0,
	          is.double(delta) && length(delta) == 1L && sigma >= 0,
	          is.integer(m) && length(m) == 1L && m >= 0L,
	          is.integer(n) && length(n) == 1L && n >= 1L,
	          is.double(init),
	          length(init) == m + n + 2L,
	          all(is.finite(init)),
	          min(init) >= 0)
	if (...length() > 0L) {
		x <- series[, 1L]
		y <- deconvolve(x = x, ...)[["value"]]
		series[, 1L] <- y[seq.int(to = length(y), length.out = length(x))]
	}
	x <- .Call(R_fastbeta, series, sigma, gamma, delta, m, n, init)
	oldClass(x) <- c("mts", "ts", "matrix", "array")
	tsp(x) <- tsp(series)
	dimnames(x) <-
		list(NULL,
		     rep(c("S", "E", "I", "R", "beta"), c(1L, m, n, 1L, 1L)))
	x
}

fastbeta.bootstrap <-
function (r,
          series, sigma = 1, gamma = 1, delta = 0,
          m = 1L, n = 1L, init, ...)
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
	          min(init) >= 0)

	## Filtering out arguments to 'seir'
	fastbeta. <-
	function (series, sigma, gamma, delta, m, n, init,
	          length.out, beta, nu, mu, stochastic, prob, delay,
	          useCompiled, ...)
	{
		## Not pretty, but neither is the (slower) alternative,
		## i.e., modifying match.call() and evaluating the result ...
		m.p <- missing(prob)
		m.d <- missing(delay)
		if (m.p && m.d)
			fastbeta(series, sigma, gamma, delta, m, n, init,
			         ...)
		else if (m.p)
			fastbeta(series, sigma, gamma, delta, m, n, init,
			         delay = delay, ...)
		else if (m.d)
			fastbeta(series, sigma, gamma, delta, m, n, init,
			         prob = prob, ...)
		else {
			if (length(prob) > 1L)
				prob <- c(rep(1, length(delay) - 1L), prob)
			fastbeta(series, sigma, gamma, delta, m, n, init,
			         prob = prob, delay = delay, ...)
		}
	}

	## Filtering out arguments to 'fastbeta'
	seir. <-
	function (length.out, beta, nu, mu, sigma, gamma, delta, m, n, init,
	          start, tol, iter.max, complete, ...)
		seir(length.out, beta, nu, mu, sigma, gamma, delta, m, n, init, ...)

	p <- m + n + 2L

	beta. <- fastbeta.(series = series,
	                   sigma = sigma, gamma = gamma, delta = delta,
	                   m = m, n = n, init = init, ...)[, p + 1L]
	nu. <- series[, 2L] # FIXME? see below
	mu. <- series[, 3L]

	length.out <- nrow(series)
	s <- as.double(seq.int(from = 0, length.out = length.out))
	beta <- approxfun(s, beta., method =   "linear", rule = 2L, ties = "ordered")
	nu   <- approxfun(s,   nu., method = "constant", rule = 2L, ties = "ordered")
	mu   <- approxfun(s,   mu., method =   "linear", rule = 2L, ties = "ordered")

	R <- simplify2array(c(list(beta.), replicate(r, simplify = FALSE, {
		X <- seir.(length.out = length.out,
		           beta = beta, nu = nu, mu = mu,
		           sigma = sigma, gamma = gamma, delta = delta,
		           m = m, n = n, init = init, ...)
		j <- p + if (ncol(X) - p == 2L) 1L:2L else 3L:2L
		series[, 1L:2L] <<- X[, j]
		fastbeta.(series = series,
		          sigma = sigma, gamma = gamma, delta = delta,
		          m = m, n = n, init = init, ...)[, p + 1L]
	})))
	oldClass(R) <- c("fastbeta.bootstrap", "mts", "ts", "matrix", "array")
	tsp(R) <- tsp(series)
	R
}

plot.fastbeta.bootstrap <-
function (x, y, level = NULL,
          col = c("#FF0000FF", "#7F7F7F40"), lwd = c(2, 1), ...)
{
	dev.hold()
	if (is.null(level))
		plot.ts(x, plot.type = "single", col = 0, lwd = 0, ...)
	else {
		stopifnot(is.double(level), length(level) == 1L, level >= 0, level <= 1)
		alpha <- 1 - level
		y <- t.default(apply(x[, -1L], 1L, quantile,
		                     probs = c(0.5 * alpha, 1 - 0.5 * alpha),
		                     names = FALSE, na.rm = TRUE))
		oldClass(y) <- oldClass(x)
		tsp(y) <- tsp(x)
		plot.ts(y, plot.type = "single", col = 0, lwd = 0, ...)
	}
	col <- rep_len(if (is.null(col)) par("col") else col, 2L)
	lwd <- rep_len(if (is.null(lwd)) par("lwd") else lwd, 2L)
	s <- as.vector(time(x))
	if (is.null(level))
		for (i in seq.int(from = 2L, length.out = ncol(x) - 1L))
			lines(s, x[, i], col = col[2L], lwd = lwd[2L])
	else {
		doPolygon <- TRUE
		n <- length(t <- s)
		if (anyNA(y)) {
			y.na <- { tmp <- is.na(y); tmp[, 1L] | tmp[, 2L] }
			if (all(y.na)) {
				warning("suppressing polygon due to NA vertices")
				doPolygon <- FALSE
			}
			else {
				i1 <-           which.min(y.na      )
				i2 <- n + 1L -  which.min(y.na[n:1L])
				i <- i1:i2
				if (any(y.na[i])) {
					warning("suppressing polygon due to NA vertices")
					doPolygon <- FALSE
				}
				else {
					n <- length(t <- t[i])
					y <- y[i, , drop = FALSE]
				}
			}
		}
		if (doPolygon)
			polygon(t[c(1L:n, n:1L)], y[c(1L:n, (n+n):(n+1L))], col = col[2L])
	}
	lines(s, x[, 1L], col = col[1L], lwd = lwd[1L])
	dev.flush()
	invisible(NULL)
}

if (FALSE) {
	nu. <-
		if (TRUE)
			## Simple but really not satisfying ...
			series[, 2L]
		else if (FALSE) {
			## Least squares _positive_ solution ... ??
			L <- suppressWarnings(
				matrix(c(0.5, 0.5, double(n)), n + 1L, n + 1L))
			L[1L] <- 1
			b <- as.double(series[, 2L])
			fn <- function (par) {
				e <- as.double(L %*% exp(par) - b)
				sum(e * e)
			}
			gr <- function (par) {
				e <- as.double(L %*% exp(par) - b)
				2 * exp(par) * colSums(e * a)
			}
			exp(optim(log(b), fn, gr, method = "BFGS")[["par"]])
		}
		else if (FALSE)
			## Can be negative, unfortunately
			suppressWarnings(
				c(1, -1) * cumsum(c(1, rep_len(c(-2, 2), n)) * series[, 2L]))
		else if (FALSE) {
			## Equivalent but slower
			L <- suppressWarnings(
				matrix(c(0.5, 0.5, double(n)), n + 1L, n + 1L))
			L[1L] <- 1
			forwardsolve(L, series[, 2L])
		}
}

fastbeta.matrix <-
function (pos,
          series, sigma = 1, gamma = 1, delta = 0,
          m = 1L, n = 1L)
{
	s <- pos
	t <- pos + 1L
	p <- m + n + 2L

	sigma <- 0.5 * m * sigma
	gamma <- 0.5 * n * gamma
	delta <- 0.5 * 1 * delta
	mu.s  <- 0.5 * 1 * series[s, 3L]
	mu.t  <- 0.5 * 1 * series[t, 3L]

	q <- c(1, sigma, gamma, delta, 0)
	a <- rep(q, c(0L, m, n, 1L, 1L))
	b <- rep(q, c(1L, m, n, 1L, 0L))

	r.0 <- (1 - mu.s - a) / (1 + mu.t + a)
	r.1 <- b / (1 + mu.t + a)

	## L := [1, 0; L0, L1]
	## y := [1; E; I; R; S]
	## iterate y = L * y = L0 + L1 * [E; I; R; S]

	L <- matrix(0, p + 1L, p + 1L)
	L[1L] <- 1

	## set band(L1,  0,  0)
	L[seq.int(from = p + 3L, by = p + 2L, length.out = p)] <- r.0
	## set band(L1, -1, -1)
	L[seq.int(from =     2L, by = p + 2L, length.out = p)] <- r.1
	## set L0[1]
	L[    2L] <- L[2L] * series[t, 1L]
	## set L0[p]
	L[p + 1L] <- (series[t, 2L] - series[t, 1L]) / (1 + mu.t)

	tmp <- L[2L, ]
	for (i in 2L:p)
		L[i + 1L, ] <- tmp <- L[i + 1L, ] + r.1[i] * tmp

	nms <- rep(c("1", "E", "I", "R", "S"), c(1L, m, n, 1L, 1L))
	dimnames(L) <- list(nms, nms)
	L
}
