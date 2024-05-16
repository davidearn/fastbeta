library(    stats, pos = "package:base", verbose = FALSE)
library(grDevices, pos = "package:base", verbose = FALSE)
library(    utils, pos = "package:base", verbose = FALSE)

library(fastbeta)
options(warn = 2L, error = if (interactive()) recover)

data(seir.E02, package = "fastbeta")
a <- attributes(seir.E02)
m <- a[["m"]]
n <- a[["n"]]

series <- cbind(seir.E02[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu")

set.seed(0L)
args <- c(list(series = series),
          a[c("sigma", "gamma", "delta", "init", "m", "n")],
          list(start = 8, end = 216))
args[["init"]] <- args[["init"]] * rlnorm(length(args[["init"]]), 0, 0.1)

L0 <- do.call(ptpi, `[[<-`(args, "complete", FALSE))
L1 <- do.call(ptpi, `[[<-`(args, "complete", TRUE ))
str(L1)

stopifnot(exprs = {
	is.list(L1)
	length(L1) == 4L
	identical(names(L1), c("value", "diff", "iter", "x"))
	is.list(L0)
	length(L0) == length(L1)
	identical(names(L0), names(L1))
	identical(L0[-4L], L1[-4L])
	is.null(L0[[4L]])
})

value <- L1[["value"]]
stopifnot(exprs = {
	is.double(value)
	length(value) == m + n + 2L
	identical(names(value), rep.int(c("S", "E", "I", "R"), c(1L, m, n, 1L)))
	!anyNA(value)
	min(value) >= 0
})

diff <- L1[["diff"]]
stopifnot(exprs = {
	is.double(diff)
	length(diff) == 1L
	!is.na(diff)
	diff >= 0
})

iter <- L1[["iter"]]
stopifnot(exprs = {
	is.integer(iter)
	length(iter) == 1L
	!is.na(iter)
	iter >= 1L
})

x <- L1[["x"]]
stopifnot(exprs = {
	is.double(x)
	is.ts(x) && inherits(x, "mts") # is.mts checks is.matrix
	identical(dim(x), c(as.integer(216 - 8 + 1), m + n + 2L, iter))
	identical(dimnames(x), list(NULL, names(value), NULL))
	!anyNA(x)
	min(x) >= 0
})

if (dev.interactive(TRUE))
	plot(x[, "S", ], plot.type = "single")

proc.time()
