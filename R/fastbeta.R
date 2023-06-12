fastbeta <-
function(Z, B, mu, gamma, S0, I0) {
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
function(Z, B, mu,
         start, a = 1L, b = length(cases),
         tol = 1e-06, iter.max = 20L, complete = FALSE) {
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
