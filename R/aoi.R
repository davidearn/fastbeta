sir.aoi <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = 1, n = max(length(R0), length(ell)),
          init = c(1 - init.infected, init.infected),
          init.infected = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, n - 1L)),
          root = NULL, aggregate = FALSE, ...)
{
	tau <- seq.int(from = from, to = to, by = by)
	stopifnot(requireNamespace("deSolve"),
	          length(tau) >= 2L, tau[1L] < tau[2L],
	          is.integer(n), length(n) == 1L, n >= 1L, n < 4096L,
	          is.double(R0), length(R0) == 1L || length(R0) == n,
	          all(is.finite(R0)), min(R0) >= 0, max(R0) > 0,
	          is.double(ell), length(ell) == 1L || length(ell) == n,
	          all(is.finite(ell)), min(ell) > 0,
	          is.double(init), length(init) == 2L,
	          all(is.finite(init)), min(init) > 0, sum(init) <= 1,
	          is.double(weights), length(weights) == n,
	          all(is.finite(weights)), min(weights) >= 0, sum(weights) > 0)
	if (!is.null(root))
	stopifnot(is.function(root),
	          !is.primitive(root),
	          all(names(formals(root)) %in% c("tau", "S", "I", "Y", "dS", "dI", "dY", "R0", "ell")))

	if (length(R0) != n)
		R0 <- rep_len(R0/n, n)
	if (length(ell) != n)
		ell <- rep_len(1, n)
	ell <- ell/sum(ell[R0 > 0])

	p <- log(init[1L])
	q <- log(init[2L])
	init <- c(p, log(weights) - log(sum(weights)) + q, q,
	          use.names = FALSE)
	if (min(init) == -Inf)
		init[init == -Inf] <- log(0x1p-64) + q

	i.S <- 1L
	i.I <- (1L + 1L):(1L + n)
	i.Y <- d <- 1L + n + 1L

	j.1 <- seq.int(from = 1L, length.out = n - 1L)
	j.2 <- seq.int(from = 2L, length.out = n - 1L)

	## D[i, j] = d(rate of change in state i)/d(state j)
	##
	## nonzero pattern in n = 4 case :
	##
	##       S I I I I Y
	##     S   3 3 3 3 .
	##     I 4 4 4 4 4 .
	##     I . 1 2 . . .
	##     I . . 1 2 . .
	##     I . . . 1 2 .
	##     Y 5 5 5 5 5 5
	##
	## where log(.) is suppressed only for pretty printing

	D <- matrix(0, d, d)
	k.1 <- seq.int(from = d + i.I[1L] + 1L, by = d + 1L,
	               length.out = n - 1L)
	k.2 <- k.1 + d
	k.3 <- i.S + (i.I - 1L) * d
	k.4 <- i.I[1L] + (c(i.S, i.I, i.I[1L]) - 1L) * d
	k.5 <- i.Y + (c(i.S, i.I, i.Y) - 1L) * d

	a.1 <- 1/ell[j.1]
	a.2 <- 1/ell[j.2]
	a.3 <- 1/ell[1L]
	a.4 <- R0/ell
	a.5 <- 1/sum(R0)
	a.6 <- rep(c(0, a.5), c(1L, n + 1L))

	gg <-
	function (t, x, theta)
	{
		log.S <- x[i.S]
		log.I <- x[i.I]
		log.Y <- x[i.Y]
		u <- sum(a.4 * exp(log.I                      ))
		v <- sum(a.4 * exp(log.I - (log.I[1L] - log.S)))
		w <- sum(a.4 * exp(log.I -  log.Y             ))
		list(c(-u,
		       v - a.3,
		       a.1 * exp(log.I[j.1] - log.I[j.2]) - a.2,
		       w * (exp(log.S) - a.5)))
	}
	Dg <-
	function (t, x, theta)
	{
		log.S <- x[i.S]
		log.I <- x[i.I]
		log.Y <- x[i.Y]
		u <- sum(uu <- a.4 * exp(log.I                      ))
		v <- sum(vv <- a.4 * exp(log.I - (log.I[1L] - log.S)))
		w <- sum(ww <- a.4 * exp(log.I -  log.Y             ))
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(log.I[j.1] - log.I[j.2]))
		D[k.3] <<- -uu
		D[k.4] <<- c(v, vv, -v + vv[1L])
		D[k.5] <<- c(w, ww, -w) * (exp(log.S) - a.6)
		D
	}
	Rg <-
	function (t, x, theta)
	{
		delayedAssign("S", exp(x[i.S]))
		delayedAssign("I", exp(x[i.I]))
		delayedAssign("Y", exp(x[i.Y]))
		root(tau = t,
		     S = S, I = I, Y = Y,
		     dS = -sum(a.4 * I) * S,
		     dI = c(sum(a.4 * I) * S - a.3 * I[1L],
		            a.1 * I[j.1] - a.2 * I[j.2]),
		     dY = sum(a.4 * I) * (S - a.5),
		     R0 = R0, ell = ell)
	}
	if (!is.null(root)) {
		call <- body(Rg)[[5L]]
		body(Rg)[[5L]] <- call[c(1L, match(names(formals(root)), names(call), 0L))]
	}

	x <- deSolve::lsoda(
		y        = init,
		times    = tau,
		func     = gg,
		parms    = NULL,
		jacfunc  = Dg,
		jactype  = "fullusr",
		rootfunc = if (!is.null(root)) Rg,
		hmax     = by,
		ynames   = FALSE,
		initfunc = NULL,
		initpar  = NULL,
		...)
	dimnames(x) <- NULL

	## FIXME: diff(t[length(t) - 1:0]) != by if terminated
	if (attr(x, "istate")[1L] < 0L)
		warning("integration terminated due to unsuccessful solver call")

	if (is.null(root)) {
		ans <- exp(x[, -1L, drop = FALSE])
		if (!aggregate)
			dimnames(ans) <- list(NULL, rep(c("S", "I", "Y"), c(1L, n, 1L)))
		else {
			ans <- cbind(seir.aggregate(ans, 0L, n),
			             rowSums(ans[, 1L + which(R0 == 0), drop = FALSE]),
			             rowSums(ans[, 1L + which(R0 >  0), drop = FALSE]),
			             deparse.level = 0L)
			dimnames(ans) <- list(NULL, c("S", "I", "Y", "I.E", "I.I"))
		}
		tsp(ans) <- c(x[c(1L, nrow(x)), 1L], 1/by)
		oldClass(ans) <- c("sir.aoi", "mts", "ts", "matrix", "array")
	} else {
		if (is.null(attr(x, "troot")))
			return(NULL)
		ans <- c(x[nrow(x), 1L], exp(x[nrow(x), -1L]))
		S <- ans[1L + i.S]
		I <- ans[1L + i.I]
		Y <- ans[1L + i.Y]
		if (!aggregate)
			names(ans) <- rep(c("tau", "S", "I", "Y"), c(1L, 1L, n, 1L))
		else {
			ans <- c(ans[1L], S, sum(I), Y,
			         sum(I[R0 == 0]), sum(I[R0 > 0]))
			names(ans) <- c("tau", "S", "I", "Y", "I.E", "I.I")
		}
		attr(ans, "curvature") <- (-sum(a.4 * I) * S) * sum(a.4 * I) + (S - a.5) * sum(a.4 * c(sum(a.4 * I) * S - a.3 * I[1L], a.1 * I[j.1] - a.2 * I[j.2]))
	}
	ans
}

summary.sir.aoi <-
function (object, tol = 1e-6, ...)
{
	stopifnot(is.double(tol), length(tol) == 1L, !is.na(tol), tol > 0)
	ans <- c(NaN, NaN)
	nms <- colnames(object)
	p <- rowSums(object[, nms == "I", drop = FALSE])
	w <- which(p > 0)
	if (length(w) == 0L || w[1L] != 1L || (end <- w[length(w)]) < 3L)
		return(ans)
	p <- log(p[1L:end])
	ph <- p[1L:(end - 1L)]
	pt <- p[2L:(end - 0L)]
	r <- (pt - ph) * tsp(object)[3L]
	rh <- r[1L:(end - 2L)]
	rt <- r[2L:(end - 1L)]
	d <- (rt - rh)/abs(rh)
	if (any(ok <- rt > 0 & d >= -tol & d < 0)) {
		rle.ok <- rle(ok)
		ptr <- c(0L, cumsum(rle.ok$lengths))
		j <- which.max(rle.ok$lengths * rle.ok$values)
		i <- (ptr[j] + 1L):ptr[j + 1L]
		ans[1L] <- rh[i[which.max(d[i])]]
	}
	if (any(ok <- rh < 0 & d >= -tol & d < 0)) {
		rle.ok <- rle(ok)
		ptr <- c(0L, cumsum(rle.ok$lengths))
		j <- which.max(rle.ok$lengths * rle.ok$values)
		i <- (ptr[j] + 1L):ptr[j + 1L]
		ans[2L] <- rt[i[which.max(d[i])]]
	}
	ans
}
