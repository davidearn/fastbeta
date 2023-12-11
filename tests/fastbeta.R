library(grDevices, pos = "package:base", verbose = FALSE)
library(    stats, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = recover)

data(sir.e02, package = "fastbeta")
a <- attributes(sir.e02)

series <- cbind(sir.e02[, c("Z.obs", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z.obs", "B", "mu")
constants <- c(S0 = sir.e02[1L, "S"],
               I0 = sir.e02[1L, "I"],
               R0 = sir.e02[1L, "R"],
               gamma = a[["gamma"]],
               delta = a[["delta"]])

X <- fastbeta(series, constants,
              prob = a[["prob"]], delay = a[["delay"]])
str(X)

stopifnot(exprs = {
	is.double(X)
	is.mts(X)
	identical(dim(X), c(nrow(sir.e02), 4L))
	identical(dimnames(X), list(NULL, c("S", "I", "R", "beta")))
	identical(tsp(X), tsp(sir.e02))
	!anyNA(X[-nrow(X), ])
	min(X, 0, na.rm = TRUE) >= 0
})

if (dev.interactive(TRUE))
	plot(X)

proc.time()
