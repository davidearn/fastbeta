sir.canonical <-
function (from = 0, to = from + 1, by = 1,
          R0, ell, n = max(length(R0), length(ell)),
          init = c(1 - p, p), p = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, n - 1L)),
          aggregate = FALSE, root = c("none", "peak"), ...)
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
	root <- match.arg(root)

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
	i.1 <- seq.int(from = 2L, length.out = n - 1L)
	i.2 <- seq.int(from = 3L, length.out = n - 1L)

	a.0 <- 1/ell[1L]
	a.1 <- 1/ell[-n]
	a.2 <- 1/ell[-1L]
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
		s.1 <- sum(a.4 * (t <- exp(x.I        )))
		s.2 <- sum(a.4 *       exp(x.I - x[2L]) )
		s.3 <- sum(a.5 *  t                     )
		list(c(-s.1 * x.S,
		       s.2 * x.S - a.0,
		       a.1 * exp(x[i.1] - x[i.2]) - a.2,
		       a.3 * x.S * s.3 - s.3))
	}
	Dg <-
	function (t, x, theta)
	{
		x.S <- x[i.S]
		x.I <- x[i.I]
		s.1 <- sum(u.1 <- a.4 * (t <- exp(x.I        )))
		s.2 <- sum(u.2 <- a.4 *       exp(x.I - x[2L]) )
		s.3 <- sum(u.3 <- a.5 *  t                     )
		D[i.S, i.S] <<- -s.1
		D[ 2L, i.S] <<-  s.2
		D[  p, i.S] <<-  a.3 * s.3
		D[i.S, i.I] <<- -u.1 * x.S
		D[ 2L, i.I] <<-  u.2 * x.S
		D[  p, i.I] <<-  a.3 * u.3 * x.S - u.3
		D[ 2L,  2L] <<- -(s.2 - u.2[1L]) * x.S
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(x[i.1] - x[i.2]))
		D
	}
	Rg <-
	switch(root,
	"peak" =
	function (t, x, theta)
		a.3 * x[i.S] - 1
	)

	x. <- deSolve::lsoda(
		y        = init,
		times    = tau,
		func     = gg,
		parms    = NULL,
		jacfunc  = Dg,
		jactype  = "fullusr",
		rootfunc = switch(root, "peak" = Rg),
		hmax     = by,
		ynames   = FALSE,
		initfunc = NULL,
		initpar  = NULL,
		...)
	dimnames(x.) <- NULL

	## FIXME: diff(t[length(t) - 1:0]) != by if terminated
	if (attr(x., "istate")[1L] < 0L)
		warning("integration terminated due to unsuccessful solver call")

	if (root == "none") {
		ans <- x.[, -1L, drop = FALSE]
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
		tsp(ans) <- c(x.[c(1L, nrow(x.)), 1L], 1/by)
		oldClass(ans) <- c("sir.canonical", "mts", "ts", "matrix", "array")
	} else {
		if (is.null(attr(x., "troot")))
			return(NULL)
		ans <- x.[nrow(x.), ]
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
		attr(ans, "curvature") <- (a.3 * (-sum(a.4 * I) * S)) * sum(a.5 * I) + (a.3 * S - 1) * sum(a.5 * c(sum(a.4 * I) * S - a.0 * I[1L], a.1 * I[-n] - a.2 * I[-1L]))
	}
	ans
}

summary.sir.canonical <-
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

seir.canonical <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = (2 * n)/(3 * n + 1),
          m = 1L, n = 1L,
          init = c(1 - p, p), p = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, m + n - 1L)),
          root = c("none", "peak"), aggregate = FALSE, ...)
{
	tau <- seq.int(from = from, to = to, by = by)
	stopifnot(requireNamespace("deSolve"),
	          length(tau) >= 2L, tau[1L] < tau[2L],
	          is.integer(m), length(m) == 1L, !is.na(m),
	          m >= 0L && m < 4096L,
	          is.integer(n), length(n) == 1L, !is.na(n),
	          n >= 1L && n < 4096L,
	          is.double(R0), length(R0) == 1L, is.finite(R0),
	          R0 > 0,
	          {
	              if (m == 0L)
	                  ell <- 0
	              TRUE
	          },
	          is.double(ell), length(ell) == 1L, is.finite(ell),
	          m == 0L || (ell > 0 && ell < 1),
	          is.double(init),
	          length(init) == 2L,
	          all(is.finite(init)),
	          init[1L] >= 0,
	          init[2L] > 0,
	          sum(init) <= 1,
	          is.double(weights),
	          length(weights) == m + n,
	          all(is.finite(weights)),
	          min(weights) >= 0,
	          sum(weights) > 0)
	root <- match.arg(root)
	p <- init[2L]
	init <- c(init[1L], log(weights) - log(sum(weights)) + log(p), p,
	          use.names = FALSE)
	if (min(init) == -Inf)
		init[init == -Inf] <- log(0x1p-64) + log(p)

	i.S <- 1L
	i.I <- (1L + m + 1L):(1L + m + n)
	i.1 <- seq.int(from = 2L, length.out = 1L + m + n - 2L)
	i.2 <- seq.int(from = 3L, length.out = 1L + m + n - 2L)

	h <- (n + 1)/2/n
	q.1 <- 1 * m/ell
	q.2 <- h * n/(1 - ell)
	q.3 <- if (m) q.1 else q.2
	q.4 <- R0 * (q.5 <- h/(1 - ell))
	a.1 <-        rep(c(q.1, q.2), c(m, n - 1L))
	a.2 <- if (m) rep(c(q.1, q.2), c(m - 1L, n)) else a.1

	## D[i, j] = d(rate of change in state i)/d(state j)
	##
	## nonzero pattern in m = 4, n = 4 case :
	##
	##       S E E E E I I I I Y
	##     S | . . . . | | | | .
	##     E | | . . . | | | | .
	##     E . | | . . . . . . .
	##     E . . | | . . . . . .
	##     E . . . | | . . . . .
	##     I . . . . | | . . . .
	##     I . . . . . | | . . .
	##     I . . . . . . | | . .
	##     I . . . . . . . | | .
	##     Y | . . . . | | | | .
	##
	## nonzero pattern in m = 0, n = 4 case :
	##
	##       S I I I I Y
	##     S | | | | | .
	##     I | | | | | .
	##     I . | | . . .
	##     I . . | | . .
	##     I . . . | | .
	##     Y | | | | | .
	##
	## where log(.) is suppressed only for pretty printing

	p <- 1L + m + n + 1L
	D <- matrix(0, p, p)
	k.1 <- seq.int(from = p + 3L, by = p + 1L, length.out = p - 3L)
	k.2 <- k.1 + p

	gg <-
	function (t, x, theta)
	{
		x.S <- x[i.S]
		x.I <- x[i.I]
		s.1 <- sum(exp(x.I        ))
		s.2 <- sum(exp(x.I - x[2L]))
		list(c(-q.4 * s.1 * x.S,
		       q.4 * s.2 * x.S - q.3,
		       a.1 * exp(x[i.1] - x[i.2]) - a.2,
		       q.4 * s.1 * x.S - q.5 * s.1))
	}
	Dg <-
	function (t, x, theta)
	{
		x.S <- x[i.S]
		x.I <- x[i.I]
		s.1 <- sum(u.1 <- exp(x.I        ))
		s.2 <- sum(u.2 <- exp(x.I - x[2L]))
		D[i.S, i.S] <<- -q.4 * s.1
		D[ 2L, i.S] <<-  q.4 * s.2
		D[  p, i.S] <<-  q.4 * s.1
		D[i.S, i.I] <<- -q.4 * u.1 * x.S
		D[ 2L, i.I] <<-  q.4 * u.2 * x.S
		D[  p, i.I] <<-  q.4 * u.1 * x.S - q.5 * u.1
		D[ 2L,  2L] <<- -q.4 * (if (m) s.2 else s.2 - 1) * x.S
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(x[i.1] - x[i.2]))
		D
	}
	Rg <-
	switch(root,
	"peak" =
	function (t, x, theta)
		q.4 * x[i.S] - q.5
	)

	x. <- deSolve::lsoda(
		y        = init,
		times    = tau,
		func     = gg,
		parms    = NULL,
		jacfunc  = Dg,
		jactype  = "fullusr",
		rootfunc = switch(root, "peak" = Rg),
		hmax     = by,
		ynames   = FALSE,
		initfunc = NULL,
		initpar  = NULL,
		...)
	dimnames(x.) <- NULL

	## FIXME: diff(t[length(t) - 1:0]) != by if terminated
	if (attr(x., "istate")[1L] < 0L)
		warning("integration terminated due to unsuccessful solver call")

	if (root == "none") {
		ans <- x.[, -1L, drop = FALSE]
		ans[, (1L + 1L):(1L + m + n)] <- exp(ans[, (1L + 1L):(1L + m + n)])
		if (aggregate) {
			ans <- seir.aggregate(ans, m, n)
			if (m > 0L)
			m <- 1L
			n <- 1L
		}
		oldClass(ans) <- c("seir.canonical", "mts", "ts", "matrix", "array")
		tsp(ans) <- c(x.[c(1L, nrow(x.)), 1L], 1/by)
		dimnames(ans) <- list(NULL, rep(c("S", "E", "I", "Y"), c(1L, m, n, 1L)))
	} else {
		if (is.null(attr(x., "troot")))
			return(c(tau = NaN, S = NaN, E = NaN, I = NaN, Y = NaN))
		r <- x.[nrow(x.), ]
		tau <- r[1L]
		S <- r[2L]
		E <- if (m > 0L) sum(E.full <- exp(r[(2L + 1L    ):(2L + m    )])) else NaN
		I <-             sum(I.full <- exp(r[(2L + 1L + m):(2L + m + n)]))
		Y <- r[2L + m + n + 1L]
		ans <- c(tau = tau, S = S, E = E, I = I, Y = Y)
		if (m > 0L)
		attr(ans, "E.full") <- E.full
		attr(ans, "I.full") <- I.full
		attr(ans, "curvature") <- -q.4 * q.4 * S * I * I + (q.4 * S - q.5) * ((if (m > 0L) q.1 * E.full[m] else q.4 * S * I) - q.2 * I.full[n])
	}
	ans
}

summary.seir.canonical <-
function (object, tol = 1e-6, ...)
{
	stopifnot(is.double(tol), length(tol) == 1L, !is.na(tol), tol > 0)
	ans <- c(NaN, NaN)
	nms <- colnames(object)
	p <- rowSums(object[, nms == "E" | nms == "I", drop = FALSE])
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
