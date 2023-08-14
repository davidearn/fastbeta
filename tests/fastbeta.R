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

n <- 250L

set.seed(0)
X1 <- sir(n, beta, nu, mu, constants, prob = prob, delay = delay)

series. <- ts(cbind(unclass(X1[, c("Z.obs", "B")]), mu = mu(0)),
              start = tsp(X1)[1L])
constants. <- constants[1:3]

## 'deconvolve' does not like NA in first row of simulation
tools::assertError(fastbeta(series., constants., prob = prob, delay = delay))

## For now, we cope ... { but why is cbind so bad with 'ts' ?? }
X1 <- window(X1, start = length(delay))
series. <- ts(cbind(unclass(X1[, c("Z.obs", "B")]), mu = mu(0)),
              start = tsp(X1)[1L])

Y <- fastbeta(series., constants., prob = prob, delay = delay)

stopifnot(exprs = {
	is.double(Y)
	is.mts(Y)
	identical(dim(Y), dim(X1))
	identical(dimnames(X1), list(NULL, c("beta", "S", "I")))
	identical(tsp(Y), tsp(X1))
	!anyNA(Y[-nrow(Y), ])
	min(Y, na.rm = TRUE) >= 0
})

if (dev.interactive(TRUE))
	plot(Y)

unclass(proc.time())
