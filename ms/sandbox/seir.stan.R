iperm <-
function (p)
{
    p[p] <- seq_along(p)
    p
}

basis.dimension <- 8L # number of basis functions
basis.order <- 4L # number of polynomial coefficients
penalty.order <- 2L # order of difference penalty

stopifnot(basis.dimension >= basis.order,
          basis.dimension > penalty.order)

T <- 32L # number of time points
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
Z. <- QU[,  j, drop = FALSE] * rep(1/sqrt(e$values[j]), each = sum(j))

## beta is unconstrained, b ~ N(0, sigma^2/lambda)

if (FALSE) {
ss <- mgcv::s(times,
              bs = "ps",
              k = basis.dimension,
              m = c(basis.order - 2L, penalty.order))
str(ss)
sc <- mgcv::smooth.construct(object = ss,
                             data = data.frame(times),
                             knots = NULL)
str(sc)
sC <- mgcv::smoothCon(object = ss,
                      data = data.frame(times),
                      knots = NULL,
                      absorb.cons = FALSE,
                      scale.penalty = FALSE)[[1L]]
str(sC)
sR <- mgcv::smooth2random(object = sC,
                          vnames = "",
                          type = 2L)
str(sR)
}
