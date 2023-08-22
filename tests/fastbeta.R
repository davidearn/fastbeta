library(grDevices, pos = "package:base", verbose = FALSE)
library(    stats, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = recover)

data(sir.e01, package = "fastbeta")
a <- attributes(sir.e01)

series <- cbind(sir.e01[, c("Z.obs", "B")], mu = a[["mu"]](0))
constants <- c(gamma = a[["gamma"]],
               S0 = sir.e01[1L, "S"],
               I0 = sir.e01[1L, "I"])

X <- fastbeta(series, constants, prob = a[["prob"]], delay = a[["delay"]])
str(X)

stopifnot(exprs = {
	is.double(X)
	is.mts(X)
	identical(dim(X), c(nrow(sir.e01), 3L))
	identical(dimnames(X), list(NULL, c("beta", "S", "I")))
	identical(tsp(X), tsp(sir.e01))
	!anyNA(X[-nrow(X), ])
	min(X, na.rm = TRUE) >= 0
})

if (dev.interactive(TRUE))
	plot(X)

unclass(proc.time())
