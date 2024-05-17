if (!(requireNamespace("adaptivetau") && requireNamespace("deSolve")))
	q("no")

library(fastbeta)
options(warn = 2L, error = if (interactive()) utils::recover)

beta <- function (t, a = 1e-01, b = 1e-05) b * (1 + a * sinpi(t / 26))
nu   <- function (t) 1e+03
mu   <- function (t) 1e-03

sigma <- 0.5
gamma <- 0.5
delta <- 0

m <- 1L
n <- 1L
p <- m + n + 2L
init <- fastbeta:::seir.ee(beta(0), nu(0), mu(0), sigma, gamma, delta, m, n)

stopifnot(exprs = {
	is.double(init)
	length(init) == p
	!anyNA(init)
	min(init) >= 0
	all.equal(sum(init), nu(0) / mu(0))
})

init <- trunc(init)

length.out <- 250L
prob <- 0.1
delay <- diff(stats::pgamma(0L:8L, 2.5))

## At the very least, these should not signal warnings or
## errors unexpectedly, and the compiled and uncompiled
## code should generate equal results

seir2 <-
function (stochastic, useCompiled)
{
	set.seed(0L)
	seir(length.out, beta, nu, mu, sigma, gamma, delta, init, m, n,
	     stochastic, prob, delay, useCompiled)
}

L <- .mapply(seir2,
             list(stochastic  = c(FALSE,  TRUE, FALSE, TRUE),
                  useCompiled = c(FALSE, FALSE,  TRUE, TRUE)),
             NULL)
dim(L) <- c(2L, 2L)

X <- L[[1L, 1L]]

stopifnot(exprs = {
	all.equal(L[, 1L], L[, 2L], tolerance = 1/sum(init))
	is.double(X)
	stats::is.mts(X)
	identical(dim(X), c(length.out, p + 3L))
	identical(dimnames(X), list(NULL, rep.int(c("S", "E", "I", "R", "Z", "B", "Z.obs"), c(1L, m, n, 1L, 1L, 1L, 1L))))
	identical(stats::tsp(X), c(0, length.out - 1, 1))
	!anyNA(X[-1L, ])
	min(0, X, na.rm = TRUE) >= 0
})

if (grDevices::dev.interactive(TRUE))
	plot(X)

tools::assertError(seir(0L, beta, nu, mu, sigma, gamma, delta, init, m, n))

proc.time()
