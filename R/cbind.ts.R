cbind.ts.dimnames <-
function (..., deparse.level = 1)
{
	stopifnot(is.numeric(deparse.level),
	          length(deparse.level) == 1L,
	          0L <= deparse.level, deparse.level < 3L)
	deparse.level <- as.integer(deparse.level)
	args <- list(NULL, ...)
	subs <- substitute(.(...))
	tags <- names(subs)
	tags.nzchar <- if (is.null(tags)) logical(length(subs)) else nzchar(tags)
	have.colnames <- FALSE
	for (i in seq_along(args))
		if (!is.null(x <- args[[i]])) {
			if (if (length(dim(x)) == 2L)
			    	!is.null(dimnames(x)[[2L]])
			    else tags.nzchar[i] || deparse.level == 2L || (deparse.level == 1L && is.symbol(subs[[i]]))) {
				have.colnames <- TRUE
				break
			}
		}
	if (have.colnames) {
		colnames <- vector("list", length(args))
		for (i in seq_along(args))
			if (!is.null(x <- args[[i]])) {
				colnames[[i]] <-
				if (length(d <- dim(x)) == 2L) {
					if (is.null((dn <- dimnames(x))[[2L]]))
						character(d[2L])
					else dn[[2L]]
				}
				else if (tags.nzchar[i])
					tags[i]
				else if (deparse.level == 2L || (deparse.level == 1L && is.symbol(subs[[i]])))
					deparse(subs[[i]])[1L]
				else ""
			}
		list(NULL, unlist(colnames, recursive = FALSE, use.names = FALSE))
	}
}

cbind.ts <-
function (..., deparse.level = 1)
{
	X <- cbind(...)
	if (!is.mts(X))
		stop(gettextf("no time series object among arguments matching '%s'",
		              "..."),
		     domain = NA)
	dimnames(X) <- cbind.ts.dimnames(..., deparse.level = deparse.level)
	X
}
