sir.aoi <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = 1, eps = 0, n = max(length(R0), length(ell)),
          init = c(1 - init.infected, init.infected),
          init.infected = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, n - 1L)),
          F = function (x) 1, Fargs = list(),
          H = identity, Hargs = list(),
          method = c("lsoda", "lsode", "vode", "daspk", "radau"),
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
	          is.double(eps), length(eps) == 1L,
	          is.finite(eps), eps >= 0,
	          is.double(init), length(init) == 2L,
	          all(is.finite(init)), init[1L] >= 0, init[2L] > 0,
	          sum(init) <= 1,
	          is.double(weights), length(weights) == n,
	          all(is.finite(weights)), min(weights) >= 0,
	          sum(weights) > 0,
	          is.function(F),
	          !is.null(formals(F)),
	          names(formals(F))[1L] != "...",
	          is.list(Fargs),
	          is.null(names(Fargs)) || !any(names(Fargs) == names(formals(F))[1L]),
	          is.function(H),
	          !is.null(formals(H)),
	          names(formals(H))[1L] != "...",
	          is.list(Hargs),
	          is.null(names(Hargs)) || !any(names(Hargs) == names(formals(H))[1L]))
	method <- if (is.list(method) && is.object(method) && inherits(method, "rkMethod")) { method.rk <- method; "rk" } else match.arg(method)
	if (!is.null(root))
	stopifnot(method %in% c("lsoda", "lsode", "radau"),
	          is.function(root),
	          !is.null(formals(root)),
	          all(names(formals(root)) %in% c("tau", "S", "I", "Y", "dS", "dI", "dY", "R0", "ell")))

	F. <- as.function(c(alist(. =), list(as.call(c(expression(F, .), Fargs)))))
	H. <- as.function(c(alist(. =), list(as.call(c(expression(H, .), Hargs)))))
	Fprime  <- F ; body(Fprime) <- D(body(F), names(formals(F))[1L])
	Hprime  <- H ; body(Hprime) <- D(body(H), names(formals(H))[1L])
	Fprime. <- F.; body(Fprime.)[[1L]] <- quote(Fprime)
	Hprime. <- H.; body(Hprime.)[[1L]] <- quote(Hprime)

	if (length(R0) != n)
		R0 <- rep_len(R0/n, n)
	if (length(ell) != n)
		ell <- rep_len(1, n)
	ell <- ell/sum(ell[R0 > 0])

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
	##     S 3 3 3 3 3 .
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
	k.3 <- i.S     + (c(i.S, i.I         ) - 1L) * d
	k.4 <- i.I[1L] + (c(i.S, i.I, i.I[1L]) - 1L) * d
	k.5 <- i.Y     + (c(i.S, i.I, i.Y    ) - 1L) * d

	a.1 <- 1/ell[j.1]
	a.2 <- 1/ell[j.2] + eps
	a.3 <- 1/ell[1L] + eps
	a.4 <- R0/ell
	a.5 <- 1/sum(R0)
	a.6 <- eps

	q <- log(init[2L])
	init <- c(init[1L], log(weights) - log(sum(weights)) + q, q,
	          use.names = FALSE)
	if (min(init) == -Inf) {
		init[init == -Inf] <- log(0x1p-64) + q
		init[i.Y] <- log(sum(exp(init[i.I])))
	}

	gg <-
	function (t, x, theta)
	{
		    S <- x[i.S]
		log.I <- x[i.I]
		log.Y <- x[i.Y]
		f <- a.4 * F.(t); h <- H.(S)
		u <- sum(f * exp(log.I            ))
		v <- sum(f * exp(log.I - log.I[1L]))
		w <- sum(f * exp(log.I - log.Y    ))
		list(c(a.6 * (1 - S) - u * h,
		       v * h - a.3,
		       a.1 * exp(log.I[j.1] - log.I[j.2]) - a.2,
		       w * (h - a.5/F.(t))))
	}
	Dg <-
	function (t, x, theta)
	{
		    S <- x[i.S]
		log.I <- x[i.I]
		log.Y <- x[i.Y]
		f <- a.4 * F.(t); h <- H.(S); hh <- Hprime.(S)
		u <- sum(uu <- f * exp(log.I            ))
		v <- sum(vv <- f * exp(log.I - log.I[1L]))
		w <- sum(ww <- f * exp(log.I - log.Y    ))
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(log.I[j.1] - log.I[j.2]))
		D[k.3] <<- -c(u * hh + a.6, uu * h)
		D[k.4] <<-  c(v * hh, vv * h, -(v - vv[1L]) * h)
		D[k.5] <<-  c(w * hh, ww * (h - a.5), -w * (h - a.5))
		D
	}
	Rg <-
	function (t, x, theta)
	{
		delayedAssign("S",     x[i.S] )
		delayedAssign("I", exp(x[i.I]))
		delayedAssign("Y", exp(x[i.Y]))
		root(tau = t,
		     S = S, I = I, Y = Y,
		     dS = a.6 * (1 - S) - sum(a.4 * I) * F.(t) * H.(S),
		     dI = c(sum(a.4 * I) * F.(t) * H.(S) - a.3 * I[1L],
		            a.1 * I[j.1] - a.2 * I[j.2]),
		     dY = sum(a.4 * I) * F.(t) * (H.(S) - a.5),
		     R0 = R0, ell = ell)
	}
	if (!is.null(root)) {
		call <- body(Rg)[[5L]]
		body(Rg)[[5L]] <- call[c(1L, match(names(formals(root)), names(call), 0L))]
	}

	args <- c(list(y = init, times = tau, func = gg, parms = NULL),
	          switch(method, "rk" = list(method = method.rk), list(jacfunc = Dg, jactype = "fullusr")),
	          switch(method, "lsoda" =, "lsode" =, "radau" = if (!is.null(root)) list(rootfun = Rg)),
	          list(ynames = FALSE, ...))
	x <- do.call(eval(substitute(deSolve::name, list(name = as.name(method)))), args)
	dimnames(x) <- NULL

	## FIXME: x[nrow(x), 1L] - x[nrow(x) - 1L, 1L] != by if terminated
	if (attr(x, "istate")[1L] < 0L)
		warning("integration terminated due to unsuccessful solver call")

	if (is.null(root)) {
		tau <- x[, 1L]
		ans <- x[, -1L, drop = FALSE]
		S <-     ans[, i.S, drop = FALSE]
		I <- exp(ans[, i.I, drop = FALSE])
		Y <- exp(ans[, i.Y, drop = FALSE])
		if (!aggregate) {
			ans[, i.I] <- I
			ans[, i.Y] <- Y
			dimnames(ans) <- list(NULL, rep(c("S", "I", "Y"), c(1L, n, 1L)))
		} else {
			ans <- cbind(S, rowSums(I), Y,
			             rowSums(I[, R0 == 0, drop = FALSE]),
			             rowSums(I[, R0 >  0, drop = FALSE]),
			             rowSums(rep(a.4, each = nrow(ans)) * I) * F.(tau),
			             rowSums(rep(a.4, each = nrow(ans)) * I) * F.(tau) * H.(S),
			             deparse.level = 0L)
			dimnames(ans) <- list(NULL, c("S", "I", "Y", "I.E", "I.I", "foi", "inc"))
		}
		tsp(ans) <- c(tau[c(1L, length(tau))], 1/by)
		oldClass(ans) <- c("sir.aoi", "mts", "ts", "matrix", "array")
		attr(ans, "eps") <- eps
	} else {
		if (is.null(attr(x, "troot")))
			return(NULL)
		tau <- x[nrow(x), 1L]
		ans <- x[nrow(x), ]
		S <-     ans[1L + i.S]
		I <- exp(ans[1L + i.I])
		Y <- exp(ans[1L + i.Y])
		if (!aggregate) {
			ans[1L + i.I] <- I
			ans[1L + i.Y] <- Y
			names(ans) <- rep(c("tau", "S", "I", "Y"), c(1L, 1L, n, 1L))
		} else {
			ans <- c(tau, S, sum(I), Y,
			         sum(I[R0 == 0]),
			         sum(I[R0 >  0]),
			         sum(a.4 * I) * F.(tau),
			         sum(a.4 * I) * F.(tau) * H.(S))
			names(ans) <- c("tau", "S", "I", "Y", "I.E", "I.I", "foi", "inc")
		}
		attr(ans, "curvature") <- (sum(a.4 * I) * Fprime.(tau) * (H.(S) - a.5)) + (sum(a.4 * I) * F.(tau) * Hprime.(S) * (a.6 * (1 - S) - sum(a.4 * I) * H.(S))) + (sum(a.4 * c(sum(a.4 * I) * H.(S) - a.3 * I[1L], a.1 * I[j.1] - a.2 * I[j.2])) * F.(tau) * (H.(S) - a.5))
	}
	ans
}

summary.sir.aoi <-
function (object, name = "Y", tol = 1e-6, ...)
{
	stopifnot(attr(object, "eps") == 0,
	          is.character(name), length(name) == 1L,
	          name %in% colnames(object), name != "S",
	          is.double(tol), length(tol) == 1L, tol > 0)
	ans <- c(NaN, NaN)
	nms <- colnames(object)
	p <- rowSums(object[, colnames(object) == name, drop = FALSE])
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
