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

## At the very least, these should not signal warnings or
## errors unexpectedly, and the compiled and uncompiled
## code should generate equal results

X00 <- sir(n, beta, nu, mu, constants, stochastic = FALSE,
           prob, delay, useCompiled = FALSE)

X01 <- sir(n, beta, nu, mu, constants, stochastic = FALSE,
           prob, delay, useCompiled =  TRUE)

set.seed(0)
X10 <- sir(n, beta, nu, mu, constants, stochastic =  TRUE,
           prob, delay, useCompiled = FALSE)

set.seed(0)
X11 <- sir(n, beta, nu, mu, constants, stochastic =  TRUE,
           prob, delay, useCompiled =  TRUE)

stopifnot(exprs = {
	all.equal(X00, X01)
	all.equal(X10, X11)
	is.double(X11)
	is.mts(X11)
	identical(dim(X11), c(n + 1L, 6L))
	identical(dimnames(X11), list(NULL, c("S", "I", "R", "B", "Z", "Z.obs")))
	identical(tsp(X11), c(0, n, 1))
	!anyNA(X11[-1L, ])
	min(X11, na.rm = TRUE) >= 0
})

if (dev.interactive(TRUE))
	plot(X11)

tools::assertError(sir(0L, beta, nu, mu, constants))

proc.time()
