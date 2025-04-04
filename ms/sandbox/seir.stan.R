library(mgcv)
invertPerm <- function(p) { p[p] <- seq_along(p); p }

basis.dimension <- 8L
basis.order <-  2L
penalty.order <- 2L

T <- 32L
K <- basis.dimension + basis.order + 2L # for P-spline

ss <- s(time,
        bs = "ps",
        k = basis.dimension,
        m = c(basis.order, penalty.order))
sc <- smooth.construct(ss, data = data.frame(time = 0:T), knots = NULL)

S <- sc[["S"]][[1L]]
X <- sc[["X"]]
F <- qr(X)
p <- invertPerm(F$pivot)
P.Rinv <- backsolve(F$qr, diag(dim(F$qr)[2L]))[p, ]

e <- eigen(crossprod(P.Rinv, S %*% P.Rinv))
abs(e$values/e$values[1L]) >= sqrt(.Machine$double.eps)


U <- e$vectors



d <- e$values

Xpri <- qr.Q(F) %*% U


sc2 <- smoothCon(ss, data=data.frame(time = 0:10), knots=NULL,
                 absorb.cons=TRUE)[[1L]]
sr.sc2 <- smooth2random(sc2, vnames="", type=2L)

## https://doi.org/10.1198/016214504000000980
