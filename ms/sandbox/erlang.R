##  Nondimensionalized interface to fastbeta::seir also computing some
##  summary statistics
##
##  [in]  from, to, by  define a sequence of increasing time points in
##                      units of the mean generation interval.
##  [in]            R0  basic reproduction number.
##  [in]           ell  ratio of the mean latent period and mean
##                      generation interval, equal to 0 if m=0 and in
##                      (0,1) otherwise.
##  [in]          m, n  number of exposed (>=0), infectious (>=1)
##                      compartments.
##  [in]          init  initial state c(X, Y) where X=S/N, Y=(E+I)/N,
##                      0 <= X, Y, X+Y <= 1.
##  [in]            y0  initial prevalence determining 'init' if 'init'
##                      is unset.  It is an error to set both 'init' and
##                      'y0'.
##  [in]            yw  numeric vector of length m+n containing weights
##                      such that yw/sum(yw) is the initial distribution
##                      of infecteds over the m+n infected compartments.
## [out] tau, x, ye, y  define time series of X=S/N, YE=E/N, Y=(E+I)/N.
## [out]          rate  numeric vector of length 2 containing the
##                      exponential rate of change at the start and end
##                      of the time series.  The first element is
##                      positive if head(Y) is increasing, else NaN.
##                      The second element is negative if tail(Y) is
##                      decreasing, else NaN.
## [out]          peak  numeric vector of length 4 of the form
##                      c(tau, x, ye, y) specifying the time in units of
##                      the mean generation interval at which Y attains
##                      a local maximum together with the corresponding
##                      values of X, YE, and Y.  If Y is monotone, then
##                      all of the elements are NaN.  If m=0, then the
##                      third element is NaN.
## [out]     peak.info  list of length 3 of the form
##                      list(ye.sub, yi.sub, curvature) storing E[i]/N,
##                      I[j]/N and Y''=((E+I)/N)'' at time peak["tau"].
##                      If Y is monotone, then all of the elements are
##                      NULL.  If m=0, then the first element is NULL.
## [out]          call  original function call after matching and
##                      evaluation of arguments.
erlang <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = if (m == 0L) 0 else (2 * n)/(3 * n + 1),
          m = 1L, n = 1L, init = c(1 - y0, y0),
          y0 = 0x1p-64, yw = rep(c(1, 0), c(1L, m + n - 1L)),
          useCompiled = TRUE) {
    call <- match.call()
    call <- as.call(c(list(quote(erlang)), mget(names(call)[-1L])))
    ## NB: This function works in units of the mean generation interval.
    ##     fastbeta::seir works in units of the observation interval.
    tau <- seq.int(from = from, to = to, by = by)
    m <- as.integer(m)
    n <- as.integer(n)
    stopifnot(length(tau) >= 2L, tau[1L] < tau[2L],
              m >= 0L,
              n >= 1L,
              is.double(R0), length(R0) == 1L,
              R0 >= 0 && R0 < Inf,
              is.double(ell), length(ell) == 1L,
              if (m == 0L) ell == 0 else ell > 0 && ell < 1,
              missing(init) || missing(y0),
              is.double(init),
              length(init) == 2L,
              all(is.finite(init)),
              min(init) >= 0,
              sum(init) <= 1,
              is.double(yw),
              length(yw) == m + n,
              all(is.finite(yw)),
              min(yw) >= 0,
              sum(yw) > 0)
    beta <- function (t) {}
    sigma <- by/ell
    gamma <- by/(1 - ell) * (n + 1)/(2 * n)
    body(beta, envir = emptyenv()) <- R0 * gamma
    init. <- c(init[1L], yw/sum(yw) * init[2L], 1 - sum(init))
    if (min(init.[-1L]) == 0) {
        ## fastbeta::seir handles E[i], I[j], R on logarithmic scale
        warning(warningCondition("setting zero-valued E[i], I[j], R to 2^-256",
                                 class = "zeroReplacedWarning"))
        init.[-1L][init.[-1L] == 0] <- 0x1p-256 # == 2^-256
        init. <- init./sum(init.)
    }
    out <- fastbeta::seir(length.out = length(tau),
                          beta = beta, sigma = sigma, gamma = gamma,
                          m = m, n = n, init = init.,
                          stochastic = FALSE, useCompiled = useCompiled)
    length(tau) <- nrow(out) # in case of early termination
    x  <- as.double(out[, 1L])
    if (m == 0L)
    y  <- rowSums(out[, (1L + 1L):(1L +     n), drop = FALSE])
    else {
    ye <- rowSums(out[, (1L + 1L):(1L + m    ), drop = FALSE])
    y  <- rowSums(out[, (1L + 1L):(1L + m + n), drop = FALSE])
    }
    ans <- list(tau = tau, x = x, ye = if (m > 0L) ye, y = y,
                rate = NULL, peak = NULL, peak.info = NULL,
                from = from, to = to, by = by,
                R0 = R0, ell = ell, m = m, n = n, init = init, yw = yw,
                call = call)
    end <-
        if (y[length(y)] > 0)
            length(y)
        else {
            ## End at last nonzero Y as Y underflowed to zero
            ## MJ: option to not invert logarithm in fastbeta::seir (?)
            w <- which(y > 0)
            max(2L, w[length(w)])
        }
    ## Return early if head(Y) nonincreasing or tail(Y) nondecreasing
    if (y[1L] >= y[2L] || y[end - 1L] <= y[end])
        return(ans)
    w <- which.max(y)
    peak.info <- list(ye.sub = NULL, yi.sub = NULL, curvature = NULL)
    ## y''
    ## = (beta S I - n gamma I[n])'
    ## = beta (S I' + S' I) - n gamma I[n]'
    if (m == 0L) {
    peak.info[["yi.sub"]] <- i.. <- as.double(out[w, (1L +     1L):(1L +     n)])
    peak.info[["curvature"]] <-
        {
            s <- x[w]
            i <- sum(i..)
            beta <- beta(0)/by
            gamma <- gamma/by
            beta * (s * (beta * s * i - n * gamma * i..[n]) + (-beta * s * i) * i) -
                n * gamma * ((if (n == 1L) beta * s * i else n * gamma * i..[n - 1L]) - n * gamma * i..[n])
        }
    } else {
    peak.info[["ye.sub"]] <- e.. <- as.double(out[w, (1L +     1L):(1L + m    )])
    peak.info[["yi.sub"]] <- i.. <- as.double(out[w, (1L + m + 1L):(1L + m + n)])
    peak.info[["curvature"]] <-
        {
            s <- x[w]
            i <- sum(i..)
            beta <- beta(0)/by
            gamma <- gamma/by
            beta * (s * (m * sigma * e..[m] - n * gamma * i..[n]) + (-beta * s * i) * i) -
                n * gamma * ((if (n == 1L) m * sigma * e..[m] else n * gamma * i..[n - 1L]) - n * gamma * i..[n])
        }
    }
    peak <- c(tau = tau[w], x = x[w], ye = NaN, y = NaN)
    if (m == 0L)
    peak[["y" ]] <- sum(i..)
    else {
    peak[["ye"]] <- sum(e..)
    peak[["y" ]] <- sum(e.., i..)
    }
    rate <- c(NaN, (log(y[end]) - log(y[end - 1L]))/by)
    if (w >= 6L) {
        z <- as.double(out[, 1L + m + n + 2L])
        w. <- which.max(z)
        tau. <- tau[4L:w.]
        z. <- z[5L:w.]
        size <- x[1L] + Wp(-R0 * x[1L] * exp(-R0 * (x[1L] + y[1L])))/R0
        f <- function(par) sum((diff(size/(1 + exp(-exp(par) * (tau. - tau.[w. - 3L])))) - z.)^2)
        rate[1L] <- exp(optimize(f, c(-100, 100))[["minimum"]])
    }
    ans[["rate"     ]] <- rate
    ans[["peak"     ]] <- peak
    ans[["peak.info"]] <- peak.info
    ans
}

Wp <-
function(z, iter.max = 10L, eps = 100 * .Machine[["double.eps"]]) {
    stopifnot(is.double(z), all(is.finite(r <- range(0, z))),
              r[1L] >= -exp(-1), r[2L] <= 0,
              is.integer(iter.max), length(iter.max) == 1L, !is.na(iter.max),
              is.double(eps), length(eps) == 1L, is.finite(eps), eps >= 0)
    w <- sqrt(2 * exp(1) * z + 2) - 1
    iter <- 0L
    done <- FALSE
    while (iter < iter.max) {
        p <- exp(w)
        t <- w * p - z
        f <- w != -1
        w <- w - f * t / (p * (w + f) - 0.5 * (w + 2) * t / (w + f))
        aok <- all(ok <- is.finite(t) & is.finite(w)) # ever not TRUE?
        if (all(abs(if (aok) t else t[ok]) <= eps * (1 + abs(if (aok) w else w[ok])))) {
            done <- TRUE
            break
        }
        iter <- iter - 1L
    }
    if (!done)
        warning(gettextf("iteration limit (%d) reached in Lambert W approximation",
                         iter.max),
                domain = NA)
    w
}
