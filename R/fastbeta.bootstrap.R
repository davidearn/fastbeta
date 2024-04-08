fastbeta.bootstrap <-
function (r, series, constants, m = 0L, n = 1L, ...)
{
	stopifnot(is.integer(m),
	          length(m) == 1L,
	          m >= 0L,
	          is.integer(n),
	          length(n) == 1L,
	          n >= 1L,
	          is.mts(series),
	          is.double(series),
	          ncol(series) == 3L,
	          min(0, series, na.rm = TRUE) >= 0,
	          is.double(constants),
	          length(constants) == m + n + 5L,
	          all(is.finite(constants)),
	          min(constants) >= 0)

	## Filtering out arguments to 'seir'
	fastbeta. <-
	function (series, constants, m, n,
	          length.out, beta, nu, mu, stochastic, prob, delay,
	          useCompiled, ...)
	{
		## Not pretty, but neither is the (slower) alternative,
		## i.e., modifying match.call() and evaluating the result ...
		m.p <- missing(prob)
		m.d <- missing(delay)
		if (m.p && m.d)
			fastbeta(series, constants, m, n, ...)
		else if (m.p)
			fastbeta(series, constants, m, n, delay = delay, ...)
		else if (m.d)
			fastbeta(series, constants, m, n, prob = prob, ...)
		else {
			if (length(prob) > 1L)
				prob <- c(rep.int(1, length(delay) - 1L), prob)
			fastbeta(series, constants, m, n, prob = prob, delay = delay, ...)
		}
	}

	## Filtering out arguments to 'fastbeta'
	seir. <-
	function (length.out, beta, nu, mu, constants, m, n,
	          start, tol, iter.max, complete, ...)
		seir(length.out, beta, nu, mu, constants, m, n, ...)

	p <- m + n + 2L

	beta. <- fastbeta.(series = series, constants = constants,
	                   m = m, n = n, ...)[, p + 1L]
	nu. <- series[, 2L] # FIXME? see below
	mu. <- series[, 3L]

	length.out <- nrow(series)
	s <- as.double(seq.int(0, length.out = length.out))
	beta <- approxfun(s, beta., method =   "linear", rule = 2L, ties = "ordered")
	nu   <- approxfun(s,   nu., method = "constant", rule = 2L, ties = "ordered")
	mu   <- approxfun(s,   mu., method =   "linear", rule = 2L, ties = "ordered")

	R <- simplify2array(c(list(beta.), replicate(r, simplify = FALSE, {
		X <- seir.(length.out = length.out, beta = beta, nu = nu, mu = mu,
		           constants = constants, m = m, n = n, ...)
		j <- p + if (ncol(X) - p == 2L) 1L:2L else 3L:2L
		series[, 1L:2L] <<- X[, j]
		fastbeta.(series = series, constants = constants,
		          m = m, n = n, ...)[, p + 1L]
	})))
	oldClass(R) <- c("fastbeta.bootstrap", oldClass(series))
	tsp(R) <- tsp(series)
	R
}

plot.fastbeta.bootstrap <-
function (x, y, level = NULL,
          col = c("#FF0000FF", "#7F7F7F40"), lwd = c(2, 1), ...)
{
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
		for (i in seq.int(2L, length.out = ncol(x) - 1L))
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
			fn <- function(par) {
				e <- as.double(L %*% exp(par) - b)
				sum(e * e)
			}
			gr <- function(par) {
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
