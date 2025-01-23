##  Nondimensionalized interface to fastbeta::seir also computing some
##  summary statistics
##
##  [in]  from, to, by  define a sequence of increasing time points in
##                      units of the mean generation interval.
##  [in]            R0  the basic reproduction number.
##  [in]  sigma, gamma  rates in reciprocal units of the mean generation
##                      interval.  At most one can be set, as the mean
##                      generation interval
##
##                          1/sigma + 2 * n * (n + 1)/gamma
##
##                      is constrained to be equal to 1.  If neither is
##                      set by the user, then sigma=2 is set internally.
##  [in]          m, n  number of exposed (>=0), infectious (>=1)
##                      compartments.
##  [in]          init  initial state c(X, Y), 0 <= X, Y, X+Y <= 1;
##                      X, Y are the susceptible, infected proportions.
##  [in]            yi  initial prevalence.  It is an error to set both
##                      'init' and 'yi'.
##  [in]          ydis  numeric vector of length m+n containing weights
##                      such that init[2] * ydis/sum(ydis) gives the
##                      initial distribution of infecteds over the m+n
##                      infected compartments.
##  [in]        smooth  a logical indicating if peak prevalence should
##                      be determined from a cubic spline fit to the
##                      time series, so as to not constrain the peak to
##                      the set of points in the time series.
##  [in]     smooth.by  a step size defining the resolution of the grid
##                      on which the cubic spline is evaluated.  It is
##                      intended that this number is smaller than 'by'.
## [out]       t, x, y  define time series of X, Y
## [out]      exponent  c(r0, r1), the asymptotic exponential rates.
##                      The first element is a positive number if
##                      head(Y) is increasing, else NaN.  The second
##                      element is a negative number if tail(Y) is
##                      decreasing, else NaN.
## [out]          peak  c(tp, yp), the time at which prevalence Y
##                      attains a local maximum and the corresponding
##                      value.  NaN means that Y' never changed sign.
## [out]          call  the original call to 'erlang' after matching and
##                      evaluation of arguments.
erlang <-
function (from = 0, to = from + 1, by = 1,
          R0,
          sigma = 1/(1 - 1/gamma * (n + 1)/(2 * n)),
          gamma = 1/((1 - 1/sigma) * (2 * n)/(n + 1)),
          m = 1L,
          n = 1L,
          init = c(1 - yi, yi),
          yi = 0x1p-10, # == 2^-10
          ydis = rep(c(1, 0), c(1L, m + n - 1L)),
          smooth = TRUE, smooth.by = by * 1e-2, ...) {
    call <- match.call(expand.dots = FALSE)
    call <- as.call(c(list(quote(erlang)),
                      mget(names(call)[-c(1L, length(call))]),
                      list(...)))
    ## NB: This function works in units of the mean generation interval.
    ##     fastbeta::seir works in units of the observation interval.
    ##     Hence: sigma -> sigma * by, gamma -> gamma * by
    t <- seq.int(from = from, to = to, by = by)
    m <- as.integer(m)
    n <- as.integer(n)
    stopifnot(length(t) >= 2L, t[1L] < t[2L],
              is.double(R0), length(R0) == 1L, R0 >= 0,
              m >= 0L,
              n >= 1L,
              if (m == 0L)
                  missing(sigma) && missing(gamma)
              else missing(sigma) || missing(gamma),
              missing(init) || missing(yi),
              is.double(init),
              length(init) == 2L,
              all(is.finite(init)),
              min(init) >= 0,
              sum(init) <= 1,
              is.double(ydis),
              length(ydis) == m + n,
              all(is.finite(ydis)),
              min(ydis) >= 0,
              sum(ydis) > 0)
    if (missing(sigma) && missing(gamma))
        sigma <- if (m == 0L) Inf else 2 # or any number greater than 1
    beta <- function (t) {}
    body(beta, envir = emptyenv()) <- R0 * gamma * by
    init <- c(init[1L], init[2L] * (ydis/sum(ydis)), 1 - sum(init))
    if (min(init[-1L]) == 0) {
        ## fastbeta::seir handles E[i], I[j], R on logarithmic scale
        warning(warningCondition("setting zero-valued E[i], I[j], R to 2^-256",
                                 class = "zeroReplacedWarning"))
        init[-1L][init[-1L] == 0] <- 0x1p-256 # == 2^-256
        init <- init/sum(init)
    }
    out <- fastbeta::seir(length.out = length(t), beta = beta,
                          sigma = sigma * by, gamma = gamma * by,
                          m = m, n = n, init = init, stochastic = FALSE)
    ## X=S/N
    x <- as.double(out[, 1L])
    ## Y=(E+I)/N
    y <- rowSums(out[, 2L:(1L + m + n), drop = FALSE])
    ans <- list(t = t, x = x, y = y,
                exponent = c(NaN, NaN), peak = c(NaN, NaN),
                call = call)
    ## NaN => nonincreasing head(Y) or nondecreasing tail(Y)
    if (max(y) <= 0)
        return(ans)
    start <-
        if (y[1L] >= y[2L])
            1L
        else {
            ## Start once exponential growth rate of Y is nonincreasing
            ## to not be misled by transient behaviour
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
    ans[["exponent"]] <-
        c(if (y.[1L] < y.[2L]) (log(y.[2L]) - log(y.[1L]))/by else NaN,
          if (y.[3L] > y.[4L]) (log(y.[4L]) - log(y.[3L]))/by else NaN)
    if (!anyNA(ans[["exponent"]])) {
        w <- which.max(y)
        ans[["peak"]] <-
            if (!smooth)
                ## Take maximum of time series
                c(t[w], y[w])
            else {
                ## Take maximum of cubic spline values over denser grid
                ss <- stats::smooth.spline(x = t, y = y, ...)
                t. <- seq.int(from = t[w - 1L], to = t[w + 1L],
                              by = smooth.by)
                y. <- stats::predict(ss, x = t.)[["y"]]
                w. <- which.max(y.)
                c(t.[w.], y.[w.])
            }
    }
    ans
}
