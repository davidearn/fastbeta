sir.aoi <- # R0=1.125, ell=4, m=1, n=1
function (from = 0, to = from + 1, by = 1,
          R0, ell, n = max(length(R0), length(ell)),
          init = c(1 - init.infected, init.infected),
          init.infected = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, n - 1L)),
          aggregate = FALSE, root = NULL, ...)
{
	tau <- seq.int(from = from, to = to, by = by)
	stopifnot(requireNamespace("deSolve"),
	          length(tau) >= 2L, tau[1L] < tau[2L],
	          is.integer(n), length(n) == 1L, n >= 1L, n < 4096L,
	          is.double(R0), length(R0) == 1L || length(R0) == n,
	          all(is.finite(R0)), min(R0) >= 0,
	          is.double(ell), length(ell) == 1L || length(ell) == n,
	          all(is.finite(ell)), min(ell) > 0,
	          is.double(init), length(init) == 2L,
	          all(is.finite(init)), init[1L] >= 0, init[2L] > 0,
	          sum(init) <= 1,
	          is.double(weights), length(weights) == n,
	          all(is.finite(weights)), min(weights) >= 0,
	          sum(weights) > 0)
	if (!is.null(root))
	stopifnot(is.function(root),
	          !is.primitive(root),
	          all(names(formals(root)) %in% c("tau", "S", "I", "Y", "dS", "dI", "dY", "R0", "ell")))

	if (length(R0) != n)
		R0 <- rep_len(R0, n)
	if (length(ell) != n)
		ell <- rep_len(ell, n)
	ell <- ell/sum(ell[R0 > 0])

	p <- init[2L]
	init <- c(init[1L], log(weights) - log(sum(weights)) + log(p), p,
	          use.names = FALSE)
	if (min(init) == -Inf)
		init[init == -Inf] <- log(0x1p-64) + log(p)

	i.S <- 1L
	i.I <- (1L + 1L):(1L + n)
	i.Y <- 1L + n + 1L
	j.1 <- seq.int(from = 1L, length.out = n - 1L)
	j.2 <- seq.int(from = 2L, length.out = n - 1L)

	a.0 <- 1/ell[1L]
	a.1 <- 1/ell[j.1]
	a.2 <- 1/ell[j.2]
	a.3 <- sum(R0)
	a.4 <- R0/ell
	a.5 <- sign(R0)

	## D[i, j] = d(rate of change in state i)/d(state j)
	##
	## nonzero pattern in n = 4 case :
	##
	##       S I I I I Y
	##     S | | | | | .
	##     I | | | | | .
	##     I . | | . . .
	##     I . . | | . .
	##     I . . . | | .
	##     Y | ~ ~ ~ ~ .
	##
	## where log(.) is suppressed only for pretty printing

	p <- 1L + n + 1L
	D <- matrix(0, p, p)
	k.1 <- seq.int(from = p + 3L, by = p + 1L, length.out = n - 1L)
	k.2 <- k.1 + p

	gg <-
	function (t, x, theta)
	{
		x.S <- x[i.S]
		x.I <- x[i.I]
		s.1 <- sum(a.4 * (q <- exp(x.I          )))
		s.2 <- sum(a.4 *       exp(x.I - x.I[1L]) )
		s.3 <- sum(a.5 *  q                       )
		list(c(-s.1 * x.S,
		       s.2 * x.S - a.0,
		       a.1 * exp(x.I[j.1] - x.I[j.2]) - a.2,
		       a.3 * x.S * s.3 - s.3))
	}
	Dg <-
	function (t, x, theta)
	{
		x.S <- x[i.S]
		x.I <- x[i.I]
		s.1 <- sum(u.1 <- a.4 * (q <- exp(x.I          )))
		s.2 <- sum(u.2 <- a.4 *       exp(x.I - x.I[1L]) )
		s.3 <- sum(u.3 <- a.5 *  q                       )
		D[i.S, i.S] <<- -s.1
		D[ 2L, i.S] <<-  s.2
		D[  p, i.S] <<-  a.3 * s.3
		D[i.S, i.I] <<- -u.1 * x.S
		D[ 2L, i.I] <<-  u.2 * x.S
		D[  p, i.I] <<-  a.3 * u.3 * x.S - u.3
		D[ 2L,  2L] <<- -(s.2 - u.2[1L]) * x.S
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(x.I[j.1] - x.I[j.2]))
		D
	}
	Rg <-
	function (t, x, theta)
	{
		delayedAssign("S",     x[i.S] )
		delayedAssign("I", exp(x[i.I]))
		delayedAssign("Y",     x[i.Y] )
		root(tau = t,
		     S = S, I = I, Y = Y,
		     dS = -S * sum(a.4 * I),
		     dI = c(S * sum(a.4 * I) - a.0 * I[1L],
		            a.1 * I[j.1] - a.2 * I[j.2]),
		     dY = (a.3 * S - 1) * sum(a.5 * I),
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
		ans <- x[, -1L, drop = FALSE]
		ans[, (1L + 1L):(1L + n)] <- exp(ans[, (1L + 1L):(1L + n)])
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
		ans <- x[nrow(x), ]
		S <- ans[2L]
		I <- exp(ans[(2L + 1L):(2L + n)])
		ans[(2L + 1L):(2L + n)] <- I
		if (!aggregate)
			names(ans) <- rep(c("tau", "S", "I", "Y"), c(1L, 1L, n, 1L))
		else {
			ans <- c(ans[1L], S, sum(I), ans[2L + n + 1L],
			         sum(I[R0 == 0]), sum(I[R0 > 0]))
			names(ans) <- c("tau", "S", "I", "Y", "I.E", "I.I")
		}
		attr(ans, "curvature") <- (a.3 * (-sum(a.4 * I) * S)) * sum(a.5 * I) + (a.3 * S - 1) * sum(a.5 * c(sum(a.4 * I) * S - a.0 * I[1L], a.1 * I[j.1] - a.2 * I[j.2]))
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
