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
## [out]yepeak, yipeak  numeric vectors of length giving E[i]/N, I[j]/N
##                      at the peak time.
## [out]          call  original function call after matching and
##                      evaluation of arguments.
erlang <-
function (from = 0, to = from + 1, by = 1,
          R0, ell = if (m == 0L) 0 else (2 * n)/(3 * n + 1),
          m = 1L, n = 1L, init = c(1 - y0, y0),
          y0 = 0x1p-64, yw = rep(c(1, 0), c(1L, m + n - 1L))) {
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
                          stochastic = FALSE)
    x  <- as.double(out[, 1L])
    if (m == 0L)
    y  <- rowSums(out[, (1L + 1L):(1L +     n), drop = FALSE])
    else {
    ye <- rowSums(out[, (1L + 1L):(1L + m    ), drop = FALSE])
    y  <- rowSums(out[, (1L + 1L):(1L + m + n), drop = FALSE])
    }
    ans <- list(tau = tau, x = x, ye = if (m > 0L) ye, y = y,
                rate = rep(NaN, 4L), peak = rep(NaN, 4L),
				yepeak = NULL, yipeak = NULL,
                from = from, to = to, by = by,
                R0 = R0, ell = ell, m = m, n = n, init = init, yw = yw,
                call = call)
    ## NaN => nonincreasing head(Y) or nondecreasing tail(Y)
    if (max(y) <= 0)
        return(ans)
    start <-
        if (y[1L] >= y[2L])
            1L
        else {
            ## Start once exponential rate of change of Y is
            ## nonincreasing to not be misled by transient behaviour
            i <- 2L; j <- 3L; r <- y[2L]/y[1L]
            j. <- length(y)
            while (j <= j. && r < (tmp <- y[j]/y[i])) {
                i <- j
                j <- j + 1L
                r <- tmp
            }
            i - 1L
        }
    end <-
        if (y[length(y)] > 0)
            length(y)
        else {
            ## End at last nonzero Y as Y underflowed to zero
            w <- which(y > 0)
            max(2L, w[length(w)])
        }
    y. <- y[c(start + 0:1, end - 1:0)]
    ans[["rate"]] <- tmp <-
        c(if (y.[1L] < y.[2L]) (log(y.[2L]) - log(y.[1L]))/by else NaN,
          if (y.[3L] > y.[4L]) (log(y.[4L]) - log(y.[3L]))/by else NaN)
    if (!anyNA(tmp)) {
        w <- which.max(y)
        ans[["peak"]] <-
            c(tau = tau[w], x = x[w], ye = if (m > 0L) ye[w] else NaN,
              y = y[w])
        if (m == 0L)
        ans[["yipeak"]] <- as.double(out[w, (1L +     1L):(1L +     n)])
        else {
        ans[["yepeak"]] <- as.double(out[w, (1L +     1L):(1L + m    )])
        ans[["yipeak"]] <- as.double(out[w, (1L + m + 1L):(1L + m + n)])
        }
    }
    ans
}
