library(grDevices, pos = "package:base", verbose = FALSE)
library(    stats, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = recover)

data(sir.e01, package = "fastbeta")
a <- attributes(sir.e01)

series <- cbind(sir.e01[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu")
constants <- c(S0 = sir.e01[1L, "S"],
               I0 = sir.e01[1L, "I"],
               R0 = sir.e01[1L, "R"],
               gamma = a[["gamma"]],
               delta = a[["delta"]])

a <- 8; b <- 216
str(L0 <- ptpi(series, constants, a = a, b = b, complete = FALSE))
str(L1 <- ptpi(series, constants, a = a, b = b, complete =  TRUE))

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
	length(value) == 3L
	identical(names(value), c("S", "I", "R"))
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
	is.ts(X) && inherits(X, "mts") # is.mts checks is.matrix
	identical(dim(X), c(as.integer(b - a + 1), 3L, iter))
	identical(dimnames(X), list(NULL, c("S", "I", "R"), NULL))
	is.finite(X)
	X >= 0
})

if (dev.interactive(TRUE))
	plot(X[, "S", ])

proc.time()
