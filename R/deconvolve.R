deconvolve <-
function(x, start, prob = 1, delay = 1,
         tol = 1, iter.max = 20L, complete = FALSE)
{
	stopifnot(exprs = {
		is.double(x)
		length(x) >= 1L
		is.double(prob)
		any(length(prob) == c(1L, length(x)))
		min(prob) >= 0
		max(prob) <= 1
		is.double(delay)
		length(delay) >= 1L
		min(delay) >= 0
		sum(delay) >  0
	})

	n <- length(x)
	d <- length(delay) - 1L
	u <- rep.int(c(0, 1, 0), c(d, n, d))
	i2  <- (d + 1L):(d + n    )
	i23 <- (d + 1L):(d + n + d)

	delay <- delay / sum(delay)
	delay.rev <- delay[(d + 1L):1L]
	q <- filter(u, delay.rev, sides = 1)[i23]

	lambda <- start
	iter <- 0L
	stol <- tol * n + 2 * sum(x)

	if (complete) {
		repeat {
			y <- filter(lambda, delay, sides = 1)[i2]
			s <- sum(y + x / (u[i2] <- y / x))
			if (iter >= iter.max || s < stol)
				break;
			lambda <- lambda * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		list(lambda = lambda,
		     chisq = (s - 2 * sum(x)) / n,
		     iter = iter)
	} else {
		lambda. <- matrix(0, length(start), iter.max + 1L)
		s. <- double(iter.max + 1L)
		repeat {
			y <- filter(lambda, delay, sides = 1)[i2]
			s.[iter + 1L] <- s <- sum(y + x / (u[i2] <- y / x))
			lambda.[, iter + 1L] <- lambda
			if (iter >= iter.max || s < stol)
				break;
			lambda <- lambda * filter(u, delay.rev, sides = 1)[i23] / q
			iter <- iter + 1L
		}
		list(lambda = lambda.,
		     chisq = (s. - 2 * sum(x)) / n,
		     iter = iter)
	}
}
