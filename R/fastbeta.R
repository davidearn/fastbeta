fastbeta <-
function(cases, births, mu, gamma, S0, I0) {
    stopifnot(exprs = {
        is.numeric(cases)
        min(0, cases, na.rm = TRUE) >= 0
        is.numeric(births)
        length(births) == length(cases)
        min(0, births, na.rm = TRUE) >= 0
        is.double(mu)
        length(mu) == length(cases)
        min(0, mu, na.rm = TRUE) >= 0
        is.double(gamma)
        length(gamma) == 1L
        gamma >= 0
        is.double(S0)
        length(S0) == 1L
        S0 >= 0
        is.double(I0)
        length(I0) == 1L
        I0 >= 0
    })
    storage.mode(cases) <- storage.mode(births) <- "double"
    .Call(R_fastbeta, cases, births, mu, gamma, S0, I0)
}

ptpi <-
function(cases, births, mu,
         start, a = 1L, b = length(cases),
         tol = 1e-06, iter.max = 20L, complete = FALSE) {
    stopifnot(exprs = {
        is.numeric(cases)
        length(cases) >= 2L
        min(0, cases, na.rm = TRUE) >= 0
        is.numeric(births)
        length(births) == length(cases)
        min(0, births, na.rm = TRUE) >= 0
        is.double(mu)
        length(mu) == length(cases)
        min(0, mu, na.rm = TRUE) >= 0
        is.double(start)
        length(start) == 1L
        start >= 0
        is.integer(a)
        length(a) == 1L
        a >= 1L
        a < length(cases)
        is.integer(b)
        length(b) == 1L
        b > a
        b <= length(cases)
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
    storage.mode(cases) <- storage.mode(births) <- "double"
    .Call(R_ptpi, cases, births, mu, start, a, b, tol, iter.max, complete)
}
