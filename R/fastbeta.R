fastbeta <-
function (Z, B, mu, gamma, S0, I0)
{
    stopifnot(exprs = {
        is.numeric(Z)
        min(0, Z, na.rm = TRUE) >= 0
        is.numeric(B)
        length(B) == length(Z)
        min(0, B, na.rm = TRUE) >= 0
        is.double(mu)
        length(mu) == length(Z)
        min(0, mu, na.rm = TRUE) >= 0
        is.double(gamma)
        length(gamma) == 1L
        gamma >= 0
        is.numeric(S0)
        length(S0) == 1L
        S0 >= 0
        is.numeric(I0)
        length(I0) == 1L
        I0 >= 0
    })
    storage.mode(Z) <- storage.mode(B) <-
        storage.mode(S0) <- storage.mode(I0) <- "double"
    as.data.frame(.Call(R_fastbeta, Z, B, mu, gamma, S0, I0))
}

ptpi <-
function (Z, B, mu,
          start, a = 1L, b = length(Z),
          tol = 1e-06, iter.max = 20L, complete = FALSE)
{
    stopifnot(exprs = {
        is.numeric(Z)
        length(Z) >= 2L
        min(0, Z, na.rm = TRUE) >= 0
        is.numeric(B)
        length(B) == length(Z)
        min(0, B, na.rm = TRUE) >= 0
        is.double(mu)
        length(mu) == length(Z)
        min(0, mu, na.rm = TRUE) >= 0
        is.numeric(start)
        length(start) == 1L
        start >= 0
        is.integer(a)
        length(a) == 1L
        a >= 0L
        a < length(Z)
        is.integer(b)
        length(b) == 1L
        b > a
        b < length(Z)
        is.double(tol)
        length(tol) == 1L
        !is.na(tol)
        is.integer(iter.max)
        length(iter.max) == 1L
        iter.max >= 1L
        is.logical(complete)
        length(complete) == 1L
        !is.na(complete)
    })
    storage.mode(Z) <- storage.mode(B) <- storage.mode(start) <- "double"
    .Call(R_ptpi, Z, B, mu, start, a, b, tol, iter.max, complete)
}

sir <- function (n, par, beta, nu, mu, stochastic = TRUE, prob = 1, delay = 1,
                 useCompiled = TRUE)
{
    ## TODO: sanity checks on arguments

    if (stochastic) {

        init <- c(S = par[["S0"]],
                  I = par[["I0"]],
                  R = par[["N0"]] - par[["S0"]] - par[["I0"]],
                  Q = 0)
        tran <- list(c(S =  1),               # birth
                     c(S = -1, I = 1, Q = 1), # infection
                     c(I = -1, R = 1),        # removal
                     c(S = -1),               # natural mortality
                     c(I = -1),               # ""
                     c(R = -1))               # ""
        if (useCompiled) {
            .Call(R_sir_ad_initialize, beta, nu, mu, par[["gamma"]])
            rate <- function (x, theta, t) .Call(R_sir_ad_rate, t, x)
            jaco <- function (x, theta, t) .Call(R_sir_ad_jaco, t, x)
        } else {
            J <- matrix(0, 4L, 6L)
            Ji <- c(1L, 5L, 6L, 10L, 13L, 18L, 23L)
            rate <- function (x, theta, t) {
                xS <- x[[1L]]
                xI <- x[[2L]]
                xR <- x[[3L]]
                beta <- theta[[1L]](t)
                nu <- theta[[2L]](t)
                mu <- theta[[3L]](t)
                gamma <- theta[[4L]]
                c(nu,
                  beta * xS * xI,
                  if (xI > 1) gamma * xI else 0,
                  mu * xS,
                  if (xI > 1) mu * xI else 0,
                  mu * xR)
            }
            jaco <- function (x, theta, t) {
                xS <- x[[1L]]
                xI <- x[[2L]]
                xR <- x[[3L]]
                beta <- theta[[1L]](t)
                nu <- theta[[2L]](t)
                mu <- theta[[3L]](t)
                gamma <- theta[[4L]]
                J[Ji] <<- c(nu, beta * xI, beta * xS, gamma, mu, mu, mu)
                J
            }
        }

        X. <- ssa.adaptivetau(
            init.values  = init,
            transitions  = tran,
            rateFunc     = rate,
            params       = if (!useCompiled) list(beta, nu, mu, par[["gamma"]]),
            tf           = n,
            jacobianFunc = jaco,
            tl.params    = list(maxtau = 0.999, extraChecks = FALSE))
        N <- dim(X.)[1L]
        i <- N - match(0L:n, as.integer(ceiling(X.[N:1L, 1L]))) + 1L
        X <- X.[i, 2L:5L, drop = FALSE]

    } else {

        init <- c(S = par["S0"],
                  logI = log(par["I0"]),
                  R = par[["N0"]] - par[["S0"]] - par[["I0"]],
                  Q = 0)
        if (useCompiled) {
            .Call(R_sir_de_initialize, beta, nu, mu, par[["gamma"]])
            rate <- "R_sir_de_rate"
            jaco <- "R_sir_de_jaco"
        } else {
            J <- matrix(0, 4L, 4L)
            Ji <- c(1L, 2L, 4L, 5L, 7L, 8L, 11L)
            rate <- function(t, x, theta) {
                xS <- x[[1L]]
                xI <- exp(x[[2L]])
                xR <- x[[3L]]
                beta <- theta[[1L]](t)
                nu <- theta[[2L]](t)
                mu <- theta[[3L]](t)
                gamma <- theta[[4L]]
                list(c(nu - beta * xS * xI - mu * xS,
                       beta * xS - gamma - mu,
                       gamma * xI - mu * xR,
                       beta * xS * xI))
            }
            jaco <- function(t, x, theta) {
                xS <- x[[1L]]
                xI <- exp(x[[2L]])
                xR <- x[[3L]]
                beta <- theta[[1L]](t)
                nu <- theta[[2L]](t)
                mu <- theta[[3L]](t)
                gamma <- theta[[4L]]
                beta.xI <- beta * xI
                beta.xS.xI <- beta.xI * xS
                J[Ji] <<- c(-beta.xI - mu, beta, beta.xI, -beta.xS.xI,
                            gamma * xI, beta.xS.xI, -mu)
                J
            }
        }

        X. <- ode(
            y = init,
            times = 0:n,
            func = rate,
            parms = if (!useCompiled) list(beta, nu, mu, par[["gamma"]]),
            jacfunc = jaco,
            jactype = "fullusr",
            hmax = 1,
            ynames = FALSE,
            dllname = if (useCompiled) "fastbeta",
            initfunc = NULL,
            initpar = NULL)
        X.[, 4L] <- exp(X.[, 4L])
        N <- dim(X.)[1L]
        X <-
            if (N < n + 1L)
                rbind(X.[, 2L:5L, drop = FALSE],
                      matrix(NA_real_, n + 1L - N, 4L))
            else X.[, 2L:5L, drop = FALSE]

    }
    head <- 1L:n
    tail <- 2L:(n + 1L)
    X[tail, 4L] <- Xt4 <- X[tail, 4L] - X[head, 4L]
    X[  1L, 4L] <- NA_real_

    m.p <- missing(prob)
    m.d <- missing(delay)
    if (doObs <- !(m.p && m.d)) {
        X <- cbind(X, 0, deparse.level = 0L)
        Xt5 <- as.integer(Xt4)
        if (!m.p)
            Xt5 <- rbinom(n, Xt5, prob)
        if (!m.d)
            Xt5 <- tabulate(rep.int(1L:n, Xt5) +
                            sample(seq.int(to = 0L, by = -1L,
                                           length.out = length(delay)),
                                   size = sum(Xt5),
                                   replace = TRUE,
                                   prob = delay),
                            n)
        X[tail, 5L] <- Xt5
        X[  1L, 5L] <- NA_real_
    }
    dimnames(X) <- list(time = NULL,
                        variable = c("S", "I", "R", "Z", if (doObs) "Zobs"))
    X
}
