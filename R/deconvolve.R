## FIXME: not clearly robust to many zeros in 'x' and/or 'delay'

deconvolve <-
function(x, prob = 1, delay = 1,
         start = 0.125 + c(double(length(delay) - 1L), x),
         tol = 1, iter.max = 20L, complete = FALSE)
{
	stopifnot(exprs = {
		is.numeric(x)
		length(x) >= 1L
		!anyNA(x)
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
		!anyNA(start)
	})

	storage.mode(x) <- storage.mode(start) <- "double"
	n <- length(x)
	d <- length(delay) - 1L
	u <- rep.int(c(0, 1, 0), c(d, n, d))
	i2  <- (d + 1L):(d + n)
	i23 <- (d + 1L):(d + n + d)

	if ((delay.sum <- sum(delay)) != 1)
		delay <- delay / delay.sum
	delay.rev <- delay[(d + 1L):1L]
	q <- filter(u, delay.rev, sides = 1)[i23]

	dx <- start
	iter <- 0L
	stol <- tol * n + 2 * sum(x)

	if (complete) {
		dx. <- matrix(0, length(start), iter.max + 1L)
		s. <- double(iter.max + 1L)
		repeat {
			y <- filter(dx, delay, sides = 1)[i2]
			s.[iter + 1L] <- s <- sum(y + x * (u[i2] <- x / y))
			dx.[, iter + 1L] <- dx
			if (iter >= iter.max || s < stol)
				break;
			dx <- dx * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		list(dx = if (missing(prob)) dx. else dx. / prob,
		     chisq = (s. - 2 * sum(x)) / n,
		     iter = iter)
	} else {
		repeat {
			y <- filter(dx, delay, sides = 1)[i2]
			s <- sum(y + x * (u[i2] <- x / y))
			if (iter >= iter.max || s < stol)
				break;
			dx <- dx * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		list(dx = if (missing(prob)) dx else dx / prob,
		     chisq = (s - 2 * sum(x)) / n,
		     iter = iter)
	}
}
