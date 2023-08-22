library(grDevices, pos = "package:base", verbose = FALSE)
library(    stats, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = recover)

data(sir.e01, package = "fastbeta")
a <- attributes(sir.e01)

series <- cbind(sir.e01[, c("Z", "B")], mu = a[["mu"]](0))

a <- 2296; b <- 2452
str(L0 <- ptpi(series, a = a, b = b, start = 4e+04, complete = FALSE))
str(L1 <- ptpi(series, a = a, b = b, start = 4e+04, complete =  TRUE))

stopifnot(exprs = {
	is.list(L0)
	length(L0) == 4L
	identical(names(L0), c("value", "delta", "iter", "X"))
	is.list(L1)
	length(L1) == length(L0)
	identical(names(L1), names(L0))
	identical(L0[-4L], L1[-4L])
	is.null(L0[[4L]])
})

value <- L1[["value"]]
stopifnot(exprs = {
	is.double(value)
	length(value) == 1L
	is.finite(value)
	value >= 0
})

delta <- L1[["delta"]]
stopifnot(exprs = {
	is.double(delta)
	length(delta) == 1L
	is.finite(delta)
})

iter <- L1[["iter"]]
stopifnot(exprs = {
	is.integer(iter)
	length(iter) == 1L
	!is.na(iter)
	iter >= 1L
})

X <- L1[["X"]]
stopifnot(exprs = {
	is.double(X)
	is.mts(X)
	dim(X)[2L] == iter
	is.null(dimnames(X))
	is.finite(X)
	X >= 0
})

if (dev.interactive(TRUE))
	plot(X)

unclass(proc.time())
