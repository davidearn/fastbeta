sir.e01 <- local({
beta <- function (t, a = 1e-01, b = 1e-05)
	b * (1 + a * cospi(t / 26))
nu <- function (t) 1e+03
mu <- function (t) 1e-03
environment(beta) <- environment(nu) <- environment(mu) <- .GlobalEnv

gamma <- 0.5
delta <- 0
S0 <- 5e+04
I0 <- 1e+03
R0 <- 1e+06 - S0 - I0
constants <- c(gamma = gamma, delta = delta, S0 = S0, I0 = I0, R0 = R0)

n <- 52L * 55L
prob <- 0.1
delay <- diff(stats::pgamma(0L:8L, 2.5))

set.seed(0L)
X <- fastbeta::sir(n, beta, nu, mu, constants,
                   prob = prob, delay = delay, epsilon = 0.0025)
X <- stats::window(X, start = n - 52L * 5L)
stats::tsp(X)[1L:2L] <- stats::tsp(X)[1L:2L] - stats::tsp(X)[1L]

structure(X,
          beta = beta, nu = nu, mu = mu, gamma = gamma, delta = delta,
          prob = prob, delay = delay)
})
