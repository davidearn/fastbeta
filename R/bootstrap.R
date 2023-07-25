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

	## Filtering out arguments to 'sir'
	fastbeta. <- function (series, constants,
	                       n, beta, nu, mu, stochastic, prob, delay,
	                       useCompiled, ...) {
		## Not pretty, but neither is the (slower) alternative,
		## i.e., modifying match.call() and evaluating the result ...
		m.p <- missing(prob)
		m.d <- missing(delay)
		if (m.p && m.d)
			fastbeta(series, constants, ...)
		else if (m.p)
			fastbeta(series, constants, delay = delay, ...)
		else if (m.d)
			fastbeta(series, constants, prob = prob, ...)
		else {
			if (length(prob) > 1L)
				prob <- c(rep.int(1, length(delay) - 1L), prob)
			fastbeta(series, constants, prob = prob, delay = delay, ...)
		}
	}

	## Filtering out arguments to 'fastbeta'
	sir. <- function (n, beta, nu, mu, constants,
	                  x, start, tol, iter.max, complete, ...)
		sir(n, beta, nu, mu, constants, ...)

	beta. <- fastbeta.(series = series, constants = constants, ...)[, 1L]
	nu. <- series[, 2L] # FIXME? see below
	mu. <- series[, 3L]

	n <- nrow(series) - 1L
	s <- as.double(0L:n)
	beta <- approxfun(s, beta., method =   "linear", rule = 2L, ties = "ordered")
	nu   <- approxfun(s,   nu., method = "constant", rule = 2L, ties = "ordered")
	mu   <- approxfun(s,   mu., method =   "linear", rule = 2L, ties = "ordered")
	constants0 <- c(constants, 0) # as removed population size is irrelevant

	R <- simplify2array(c(list(beta.), replicate(r, simplify = FALSE, {
		X <- sir.(n = n, beta = beta, nu = nu, mu = mu,
		          constants = constants0, ...)
		series[, 1L:2L] <<- X[, c(ncol(X), 4L)]
		fastbeta.(series = series, constants = constants, ...)[, 1L]
	})))
	oldClass(R) <- oldClass(series)
	tsp(R) <- tsp(series)
	R
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
