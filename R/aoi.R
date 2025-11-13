sir.aoi <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = 1, eps = 0, n = max(length(R0), length(ell)),
          init = c(1 - init.infected, init.infected),
          init.infected = .Machine[["double.neg.eps"]],
          weights = rep(c(1, 0), c(1L, n - 1L)),
          F = function (x) 1, Fargs = list(),
          H = identity, Hargs = list(),
          root = NULL, root.max = 1L, root.break = TRUE,
          aggregate = FALSE, skip.Y = FALSE,
          rtol = 1e-9, atol = rep(c(1e-9, 0), c(1L, n + !skip.Y)), ...)
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
	if (!is.null(root))
	stopifnot(is.function(root),
	          !is.null(formals(root)),
	          all(names(formals(root)) %in% c("tau", "S", "I", if (!skip.Y) "Y", "dS", "dI", if (!skip.Y) "dY", "R0", "ell")),
	          is.integer(root.max), length(root.max) == 1L, root.max >= 1L)

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
	i.Y <- 1L + n + 1L

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

	d <- if (skip.Y) 1L + n else 1L + n + 1L
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
	function (t, y, theta)
	{
		    S <- y[i.S]
		log.I <- y[i.I]
		log.Y <- y[i.Y]
		f <- F.(t); h <- H.(S)
		u <- sum(a.4 * f * exp(log.I            ))
		v <- sum(a.4 * f * exp(log.I - log.I[1L]))
		w <- sum(a.4 * f * exp(log.I - log.Y    ))
		list(c(a.6 * (1 - S) - u * h,
		       v * h - a.3,
		       a.1 * exp(log.I[j.1] - log.I[j.2]) - a.2,
		       w * (h - a.5/f)))
	}
	Dg <-
	function (t, y, theta)
	{
		    S <- y[i.S]
		log.I <- y[i.I]
		log.Y <- y[i.Y]
		f <- F.(t); h <- H.(S); hh <- Hprime.(S)
		u <- sum(uu <- a.4 * f * exp(log.I            ))
		v <- sum(vv <- a.4 * f * exp(log.I - log.I[1L]))
		w <- sum(ww <- a.4 * f * exp(log.I - log.Y    ))
		D[k.2] <<- -(D[k.1] <<- a.1 * exp(log.I[j.1] - log.I[j.2]))
		D[k.3] <<- -c(u * hh + a.6, uu * h)
		D[k.4] <<-  c(v * hh, vv * h, -(v - vv[1L]) * h)
		D[k.5] <<-  c(w * hh, ww * (h - a.5/f), -w * (h - a.5/f))
		D
	}
	Rg <-
	function (t, y, theta)
	{
		delayedAssign("S",     y[i.S] )
		delayedAssign("I", exp(y[i.I]))
		delayedAssign("Y", exp(y[i.Y]))
		delayedAssign("f", F.(t))
		delayedAssign("h", H.(S))
		root(tau = t,
		     S = S, I = I, Y = Y,
		     dS = a.6 * (1 - S) - sum(a.4 * f * I) * h,
		     dI = c(sum(a.4 * f * I) * h - a.3 * I[1L],
		            a.1 * I[j.1] - a.2 * I[j.2]),
		     dY = sum(a.4 * f * I) * (h - a.5/f),
		     R0 = R0, ell = ell)
	}

	if (!is.null(root)) {
		b <- body(Rg)
		b[[7L]] <- b[[7L]][c(1L, match(names(formals(root)), names(b[[7L]]), 0L))]
		body(Rg) <- b
	}

	if (skip.Y) {
		init <- init[-i.Y]
		D <- D[-i.Y, -i.Y, drop = FALSE]

		b <- body(gg)
		b[[4L]] <- b[[9L]] <- b[[c(10L, 2L, 5L)]] <- NULL
		body(gg) <- b

		b <- body(Dg)
		b[[4L]] <- b[[10L]] <- b[[14L]] <- NULL
		body(Dg) <- b

		b <- body(Rg)
		b[[4L]] <- NULL
		body(Rg) <- b
	}

	common <-
	function (t, y)
	{
		S <-     y[, i.S, drop = FALSE]
		I <- exp(y[, i.I, drop = FALSE])
		if (!skip.Y)
		Y <- exp(y[, i.Y, drop = FALSE])
		if (aggregate) {
			m <- nrow(y)
			f <- F.(t); h <- H.(S); ff <- Fprime.(t); hh <- Hprime.(S)
			a.1 <- rep(a.1, each = m)
			a.2 <- rep(a.2, each = m)
			a.4 <- rep(a.4, each = m)
			y <- cbind(S, rowSums(I), if (!skip.Y) Y,
			           rowSums(I[, R0 == 0, drop = FALSE]),
			           rowSums(I[, R0 >  0, drop = FALSE]),
			           rowSums(a.4 * f * I),
			           rowSums(a.4 * f * I) * h,
			           rowSums(a.4 * ff * I) * (h - a.5/f) +
			           	rowSums(a.4 * f * I) * (hh * (a.6 * (1 - drop(S)) - rowSums(a.4 * f * I) * h) + a.5 * ff/f/f) +
			           	rowSums(a.4 * f * cbind(rowSums(a.4 * f * I) * h - a.3 * I[, 1L], a.1 * I[, j.1, drop = FALSE] - a.2 * I[, j.2, drop = FALSE], deparse.level = 0L)) * (h - a.5/f),
			           deparse.level = 0L)
			dimnames(y) <- list(NULL, c("S", "I", if (!skip.Y) "Y", "I.E", "I.I", "foi", "inc", "crv"))
		} else {
			y <- cbind(S, I, if (!skip.Y) Y, deparse.level = 0L)
			dimnames(y) <- list(NULL, rep(c("S", "I", if (!skip.Y) "Y"), c(1L, n, if (!skip.Y) 1L)))
		}
		list(tau = t, state = y)
	}

	x <- deSolve::lsoda(y = init,
	                    times = tau,
	                    func = gg,
	                    parms = NULL,
	                    rtol = rtol,
	                    atol = atol,
	                    jacfunc = Dg,
	                    jactype = "fullusr",
	                    rootfunc =
	                        if (!is.null(root))
	                        	Rg,
	                    events =
	                        if (!is.null(root) && !(root.break && root.max == 1L))
	                        	list(func = function (t, y, theta, ...) y,
	                        	     root = TRUE,
	                        	     maxroot = root.max),
	                    ynames = FALSE, ...)
	ax <- attributes(x)
	attributes(x) <- ax["dim"]

	status <- ax[["istate"]][1L]
	if (status < 0L)
		warning("integration terminated due to unsuccessful solver call")
	if (status < 0L || status == 3L) {
		last <- x[nrow(x), , drop = FALSE]
		if (!last[1L, 1L] %in% tau)
		x <- x[-nrow(x), , drop = FALSE]
	}

	if (root.break && root.max > 1L &&
	    !is.null(ax[["nroot"]]) && ax[["nroot"]] == root.max &&
	    x[nrow(x), 1L] > ax[["troot"]][root.max])
		x <- x[x[, 1L] <= ax[["troot"]][root.max], , drop = FALSE]

	cx <- common(x[, 1L], x[, -1L, drop = FALSE])

	ans <- cx[[2L]]
	tsp(ans) <- c(cx[[1L]][c(1L, length(cx[[1L]]))], 1/by)
	oldClass(ans) <- c("sir.aoi", "mts", "ts", "matrix", "array")
	if (!is.null(root)) {
		attr(ans, "root.info") <-
		if (status == 3L)
			common(last[, 1L], last[, -1L, drop = FALSE])
		else if (is.null(ax[["nroot"]]))
			common(double(0L), matrix(0, 0L, d))
		else
			common(ax[["troot"]], t(ax[["valroot"]]))
	}
	attr(ans, "eps") <- eps
	ans
}

summary.sir.aoi <-
function (object, name = "Y", tol = 1e-6, ...)
{
	stopifnot(attr(object, "eps") == 0,
	          is.character(name), length(name) == 1L,
	          name %in% colnames(object),
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
