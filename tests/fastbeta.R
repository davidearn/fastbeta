library(    stats, pos = "package:base", verbose = FALSE)
library(grDevices, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = if (interactive()) recover)

data(seir.E02, package = "fastbeta")
a <- attributes(seir.E02)
m <- a[["m"]]
n <- a[["n"]]

series <- cbind(seir.E02[, c("Z.obs", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z.obs", "B", "mu")

X <- do.call(fastbeta,
             c(list(series = series),
               a[c("sigma", "gamma", "delta", "init", "m", "n", "prob", "delay")]))
str(X)

stopifnot(exprs = {
	is.double(X)
	is.mts(X)
	identical(dim(X), c(nrow(seir.E02), m + n + 3L))
	identical(dimnames(X), list(NULL, rep.int(c("S", "E", "I", "R", "beta"), c(1L, m, n, 1L, 1L))))
	identical(tsp(X), tsp(seir.E02))
	!anyNA(X[-length(X)])
	min(0, X, na.rm = TRUE) >= 0
})

if (dev.interactive(TRUE))
	plot(X)

proc.time()
