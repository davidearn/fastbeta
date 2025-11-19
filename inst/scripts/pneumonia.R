convolveGamma <-
function (x, shape, scale, tol = 1e-6, max.terms = 1024L,
          prec = 1024L, verbose = getOption("verbose")) {
    ## https://doi.org/10.1007/BF02481123
    loinga <- flint::arb_hypgeom_gamma_lower # [lo]wer [in]complete [ga]mma
    K <- as.integer(max.terms)
    if (K < 1L)
        return(double(length(x)))
    n <- max(length(shape), length(scale))
    alpha <- rep_len(shape, n)
    beta  <- rep_len(scale, n)
    o <- order(beta)
    alpha <- alpha[o]
    beta  <- beta [o]
    x <- x/beta[1L]
    C <- prod((beta[1L]/beta)^alpha)
    gamma. <- colSums(matrix(alpha * (1 - beta[1L]/beta)^rep(seq_len(K - 1L), each = n), nrow = n))
    rho <- sum(alpha)
    delta <- double(K)
    ans <- (ans. <- double(length(x))) +
        (delta[1L] <- 1) *
        as.double(loinga(rho, x, 1L, prec))
    k <- 1L
    while (k < K && any(abs(ans - ans.) > tol * abs(ans.))) {
        ans <- (ans. <- ans) +
            (delta[k + 1L] <- sum(gamma.[1L:k] * delta[k:1L])/k) *
            as.double(loinga(rho + k, x, 1L, prec))
        k <- k + 1L
        if (verbose)
            cat(sprintf("delta %.4e with term %d\n",
                        max(abs((ans - ans.)/ans.)[ans. > 0]), k))
    }
    `attr<-`(`attr<-`(C * ans, "terms", k),
             "delta", max(abs((ans - ans.)/ans.)[ans. > 0]))
}

optimShapeScale <-
function (F, L) {
    n <- diff(L$y)
    fn <-
    function (par) {
        p <- diff(F(L$x, shape = exp(par[[1L]]), scale = exp(par[[2L]])))
        -sum((n * log(p))[p > 0])
    }
    M <- optim(c(0, 0), fn)
    xlab <- "x"
    ylab <- sprintf("%s(%s, shape = %g, scale = %g)\n",
                    deparse(substitute(F)),
                    xlab,
                    exp(M$par[[1L]]),
                    exp(M$par[[2L]]))
    `class<-`(list(F = F, L = L, M = M, xlab = xlab, ylab = ylab),
              "oss")
}

.S3method("plot", "oss",
          function (x, y, ...) {
              s <- seq(from = min(x$L$x), to = max(x$L$x),
                       length.out = 256L)
              t <- x$F(s,
                       shape = exp(x$M$par[[1L]]),
                       scale = exp(x$M$par[[2L]]))
              plot(s, t, type = "l", xlab = x$xlab, ylab = x$ylab, ...)
              points(x$L$x, x$L$y/max(x$L$y))
              invisible(NULL)
          })

.S3method("print", "oss",
          function (x, ...) {
              cat(sprintf("%s\n", x$ylab))
              invisible(x)
          })


## Time from infection to symptom onset
## https://doi.org/10.1093/oxfordjournals.aje.a112781
csv <- read.csv("1979MoseBend-Figure1-Bars.csv",
                header = FALSE, col.names = c("x", "y"))
## y(x) := number of cases with delay <= x
(L1  <- data.frame(x = c(0, csv$x), y = cumsum(c(0, csv$y))))
plot(print(G1 <- optimShapeScale(pgamma, L1)))

## Time from symptom onset to death
## https://doi.org/10.1001/jama.1918.02600500012003
## Digitized by WebPlotDigitizer
csv <- read.csv("1918KeatCush-Chart2-SolidLine.csv",
                header = FALSE, col.names = c("x", "y"))
L2 <- approx(csv$x, csv$y, 0:max(csv$x))
L2$x <-            round(c(0, 1 + L2$x))
L2$y <- as.integer(round(c(0,     L2$y)))
L2$y[is.na(L2$y) | L2$x <= 2L | L2$x > 28L] <- 0L
L2$y <- cumsum(L2$y)
## y(x) := number of cases with delay <= x
(L2 <- as.data.frame(L2))
plot(print(G2 <- optimShapeScale(pgamma, L2)))

## Time from infection to death
## https://doi.org/10.1073/pnas.090295810
## Digitized by WebPlotDigitizer
## MJ: not reproducible by other means after several tries ...
csv <- read.csv("2009GoldDush-SupplementFigure1-Points.csv",
                header = FALSE, col.names = c("x", "y"))
L3 <- approx(csv$x, csv$y, 0:max(csv$x))
L3$y[is.na(L3$y) | L3$y < 0] <- 0
L3$x <- round(L3$x, digits = 0L)
L3$y <- round(L3$y, digits = 4L)
## y(x) := probability of delay == x
(L3 <- as.data.frame(L3))

d <- 255L; tt <- 0L:d; tttt <- 0L:(d + d)

cG <- convolveGamma(tt,
                    shape = exp(c(G1$M$par[[1L]],
                                  G2$M$par[[1L]])),
                    scale = exp(c(G1$M$par[[2L]],
                                  G2$M$par[[2L]])),
                    verbose = TRUE)
pG1G2 <- zapsmall(c(0, diff(cG)))

pW1 <- c(0, diff(pweibull(tt - 0.5, shape = 2.21 + 1, scale = 1.10)))
pN2 <- c(0, diff(L2$y/max(L2$y)), double(1L + d - nrow(L2)))
pW1N2 <- zapsmall(convolve(pW1, rev(pN2), type = "open"))

plot(0, 0, type = "n",
     xlim = c(0, 64), ylim = c(0, 0.12),
     xlab = "days", ylab = "probability")
lines(  tt, pG1G2, col = 1, lwd = 6)
lines(tttt, pW1N2, col = 2, lwd = 6)
lines(L3$x,  L3$y, col = 3, lwd = 6)

series <- read.csv("pneumonia.csv",
                   header = FALSE, col.names = c("date", "deaths"))
series$date <- `storage.mode<-`(as.Date(series$date), "integer")
series

d <- max(which(pG1G2 > 0), nrow(L3)) - 1L
delay <- data.frame(nday = .difftime(0:d, "days"),
                    goldstein = c(L3$y, double(1L + d - nrow(L3))),
                    gpg = pG1G2[seq_len(1L + d)])
delay

pneumonia <- list(series = series, delay = delay)
save(pneumonia, file = "pneumonia.rda", version = 3L, compress = "xz")
