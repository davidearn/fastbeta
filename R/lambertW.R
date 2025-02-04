lambertW <-
function (z, iter.max = 10L, eps = 100 * .Machine[["double.eps"]])
{
	stopifnot(is.double(z), min(z) >= -exp(1), max(z) <= 0,
	          is.integer(iter.max), length(iter.max) == 1L, !is.na(iter.max),
	          is.double(eps), length(eps) == 1L, is.finite(eps), eps >= 0)
	w <- sqrt(2 * exp(1) * z + 2) - 1
	iter <- 0L
	done <- FALSE
	while (iter < iter.max) {
		p <- exp(w)
		t <- w * p - z
		f <- w != -1
		w <- w - f * t / (p * (w + f) - 0.5 * (w + 2) * t / (w + f))
		aok <- all(ok <- is.finite(t) & is.finite(w)) # ever not TRUE?
		if (all(abs(if (aok) t else t[ok]) <= eps * (1 + abs(if (aok) w else w[ok])))) {
			done <- TRUE
			break
		}
		iter <- iter - 1L
	}
	if (!done)
		warning(gettextf("iteration limit (%d) reached in Lambert W approximation",
		                 iter.max),
		        domain = NA)
	w
}
