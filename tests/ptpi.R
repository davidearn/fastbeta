library(stats, pos = "package:base", verbose = FALSE)
library(fastbeta)
options(warn = 2L, error = recover)

## FIXME: use data/*.R
beta <- function (t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function (t) 1e+03
mu <- function (t) 1e-03

S0 <- 5e+04
I0 <- 1e+03
constants <- c(gamma = 0.5, S0 = S0, I0 = I0, R0 = 1e+06 - S0 - I0)

prob <- 0.1
delay <- diff(pgamma(0:8, 2.5))

n <- 2500L

set.seed(0)
X1 <- window(sir(n, beta, nu, mu, constants, prob = prob, delay = delay),
             start = 2250)

series. <- ts(cbind(unclass(X1[, c("Z", "B")]), mu = mu(0)),
              start = tsp(X1)[1L])

a <- 2296; b <- 2452
str(L0 <- ptpi(series., a = a, b = b, start = 4e+04, complete = FALSE))
str(L1 <- ptpi(series., a = a, b = b, start = 4e+04, complete =  TRUE))

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
