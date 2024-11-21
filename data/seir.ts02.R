## NB: seir.ts01 differs mainly in setting of 'stochastic'
seir.ts02 <- if (requireNamespace("adaptivetau")) local({

if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE))
	set.seed(NULL, "default", "default", "default")
.R.s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

stochastic <- TRUE

beta <- function (t, a = 1e-01, b = 1e-05) b * (1 + a * sinpi(t / 26))
nu   <- function (t) 1e+03
mu   <- function (t) 1e-03
environment(beta) <- environment(nu) <- environment(mu) <- globalenv()

sigma <- 0.5
gamma <- 0.5
delta <- 0

m <- 1L
n <- 1L
init <- fastbeta::seir.ee(beta(0), nu(0), mu(0), sigma, gamma, delta, m, n)
init <- trunc(init)

length.out <- 52L * 55L
prob <- 0.1
delay <- diff(stats::pgamma(0L:8L, 2.5))

set.seed(0L, "default", "default", "default")
X <- fastbeta::seir(length.out,
                    beta, nu, mu,
                    sigma, gamma, delta,
                    m, n, init,
                    stochastic = stochastic,
                    prob = prob, delay = delay,
                    epsilon = 0.002)
X <- stats::window(X, start = stats::tsp(X)[2L] - 52L * 5L / stats::tsp(X)[3L])
stats::tsp(X)[1L:2L] <- stats::tsp(X)[1L:2L] - stats::tsp(X)[1L]

assign(".Random.seed", .R.s, envir = globalenv())

structure(X,
          beta = beta, nu = nu, mu = mu,
          sigma = sigma, gamma = gamma, delta = delta,
          m = m, n = n, init = X[1L, seq_len(m + n + 2L)],
          stochastic = stochastic,
          prob = prob, delay = delay)

})
