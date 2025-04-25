## .... Set up spline stuff ............................................

iperm <- function (p) { p[p] <- seq_along(p); p }

basis.dimension <- 8L # number of basis functions
basis.order <- 4L # number of polynomial coefficients
penalty.order <- 2L # order of difference penalty

stopifnot(basis.dimension >= basis.order,
          basis.dimension > penalty.order)

T <- 250L # number of time points
times <- seq.int(from = 0, by = 1, length.out = T)

K <- basis.dimension + basis.order # number of knots
K. <- K - 2L * (basis.order - 1L) # number of internal knots
b <- 0.001 * (times[length(times)] - times[1L])
d <- 1.002 * (times[length(times)] - times[1L])/(K. - 1L)
knots <- seq.int(from = times[1L] - b - (basis.order - 1L) * d,
                 by = d, length.out = K)

X <- splines::splineDesign(knots = knots, x = times, ord = basis.order)
F <- qr(X)
S <- crossprod(diff(diag(basis.dimension), differences = penalty.order) %*%
               backsolve(F$qr, diag(dim(F$qr)[2L]))[iperm(F$pivot), ])

stopifnot(identical(dim(X), c(T, basis.dimension)),
          identical(dim(S), c(basis.dimension, basis.dimension)))

e <- eigen(S, symmetric = TRUE)
j <- rep(c(TRUE, FALSE), c(basis.dimension - penalty.order, penalty.order))
QU <- qr.Q(F) %*% e$vectors

X. <- QU[, !j, drop = FALSE]
Z. <- QU[,  j, drop = FALSE] * rep(1/sqrt(e$values[j]), each = T)


## .... Set up incidence stuff .........................................

N <- 1e+6L
m <- 1L
n <- 1L

beta.0 <- 1e-5
beta.1 <- 1e-1
nu.0 <- 1e+3
mu.0 <- nu.0/N
sigma <- 0.5
gamma <- 0.5
delta <- 0

beta <- function (t) beta.0 * (1 + beta.1 * sinpi(t / 26))
nu   <- function (t) rep(nu.0, length(t))
mu   <- function (t) rep(mu.0, length(t))

set.seed(0L)
init <- fastbeta::seir.ee(beta.0, nu.0, mu.0, sigma, gamma, delta, m, n)
series <- fastbeta::seir(T + 1L, beta, nu, mu, sigma, gamma, delta, m, n,
                         init, epsilon = 0.002)[-1L, ]


## .... Set up data ....................................................

data <- list(T = T,
             incidence = as.integer(series[, 1L + m + n + 1L + 1L]),
             birth = as.integer(series[, 1L + m + n + 1L + 2L]),
             death = mu(times),
             sigma = sigma,
             gamma = gamma,
             delta = delta,
             m = m,
             n = n,
             init = as.double(series[1L, 1L:(1L + m + n + 1L)]),
             R0 = ncol(X.),
             R1 = ncol(Z.),
             X0 = X.,
             X1 = Z.)
