deconvolve <-
function(x, prob = 1, delay = 1,
         start = c(rep.int(1, length(delay) - 1L), x),
         tol = 1, iter.max = 20L, complete = FALSE)
{
	stopifnot(exprs = {
		is.numeric(x)
		length(x) >= 1L
		min(x) >= 0
		sum(x) >  0
		is.double(delay)
		length(delay) >= 1L
		min(delay) >= 0
		sum(delay) >  0
		is.double(prob)
		any(length(prob) == c(1L, length(x) + length(delay) - 1L))
		min(prob) >= 0
		max(prob) <= 1
		is.numeric(start)
		length(start) == length(x) + length(delay) - 1L
		min(start) >= 0
		sum(start) >  0
	})

	storage.mode(x) <- storage.mode(start) <- "double"

	n <- n. <- length(x)
	d <- length(delay) - 1L

	if (min(x) == 0) {
		i <- range(which(x > 0))
		x <- x[i[1L]:i[2L]]
		i <- i + c(0L, d)
		start <- start[i[1L]:i[2L]]
		n <- length(x)
	}

	if (min(prob) == 0)
		prob[prob == 0] <- NaN

	if ((delay.sum <- sum(delay)) != 1)
		delay <- delay / delay.sum
	delay.rev <- delay[(d + 1L):1L]

	u <- rep.int(c(0, 1, 0), c(d, n, d))
	i2  <- (d + 1L):(d + n)
	i23 <- (d + 1L):(d + n + d)

	q <- filter(u, delay.rev, sides = 1)[i23]
	q0 <- if (min(q) == 0) which(q == 0) else integer(0L)

	r <- start
	iter <- 0L
	stol <- tol * n + 2 * sum(x)

	if (complete) {
		r. <- matrix(0, d + n, iter.max + 1L)
		s. <- double(iter.max + 1L)
		repeat {
			if (length(q0))
				r[q0] <- 0
			y <- filter(r, delay, sides = 1)[i2]
			r.[, iter + 1L] <- r
			s.[iter + 1L] <- s <- sum(y + x * (u[i2] <- x / y))
			if (iter >= iter.max || s < stol)
				break
			r <- r * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		if (iter < iter.max) {
			head <- seq_len(iter + 1L)
			r. <- r.[, head, drop = FALSE]
			s. <- s.[head]
		}
		if (length(q0))
			r.[q0, ] <- NaN
		if (n < n.)
			r. <- `[<-`(matrix(NaN, d + n., iter + 1L), i[1L]:i[2L], , r.)
		if (!missing(prob))
			r. <- r. / prob
		list(dx = r., chisq = (s. - 2 * sum(x)) / n, iter = iter)
	} else {
		repeat {
			if (length(q0))
				r[q0] <- 0
			y <- filter(r, delay, sides = 1)[i2]
			s <- sum(y + x * (u[i2] <- x / y))
			if (iter >= iter.max || s < stol)
				break
			r <- r * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		if (length(q0))
			r[q0] <- NaN
		if (n < n.)
			r <- `[<-`(rep.int(NaN, d + n.), i[1L]:i[2L], r)
		if (!missing(prob))
			r <- r / prob
		list(dx = r, chisq = (s - 2 * sum(x)) / n, iter = iter)
	}
}
