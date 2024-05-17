library(fastbeta)
data(sir.e01, package = "fastbeta")

set.seed(0L)
n <- nrow(sir.e01)
x <- rbinom(n, sir.e01[, "Z"], attr(sir.e01, "prob"))

nu <- attr(sir.e01, "nu")
mu <- attr(sir.e01, "mu")
constants <- c(gamma = attr(sir.e01, "gamma"),
               sir.e01[1L, c("S", "I", "R")])
prob <- attr(sir.e01, "prob")

p1 <- c(log(0.1), log(1e-05))
f1 <- function(par) {
	beta <- function (t, a = exp(par[[1L]]), b = exp(par[[2L]]))
		b * (1 + a * cospi(t/26))
	size <- sir(n, beta, nu, mu, constants, stochastic = FALSE)[, "Z"]
	-sum(dbinom(x, floor(size), prob, log = TRUE),
	     na.rm = TRUE)
}
o1 <- optim(p1, f1)

series <- cbind(x, sir.e01[, "B"], mu(0))
beta. <- fastbeta(series, constants[1L:3L], prob = prob)[, "beta"]
time. <- time(beta.)

p2 <- p1
f2 <- function(par) {
	beta <- function (t, a = exp(par[[1L]]), b = exp(par[[2L]]))
		b * (1 + a * cospi(t/26))
	e <- beta(time.) - beta.
	sum(e * e, na.rm = TRUE)
}
o2 <- optim(p2, f2)

p3 <- o2[["par"]]
f3 <- f1
o3 <- optim(p3, f3)

m <- 41L
G <- vapply(p1, function(x) x + seq.int(-2, 2, length.out = m), double(m))
R <- matrix(list(NULL), m, m)
for (i in seq_len(m))
for (j in seq_len(m))
R[i, j] <- optim(c(G[i, 1L], G[j, 2L]), f1)
