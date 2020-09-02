#' \loadmathjax
#' Define a seasonally forcing function for simulations
#'
#' @description
#' Defines a seasonal forcing function with environmental noise
#' that can be assigned to the argument `beta` of [make_data()].
#'
#' @details
#' A noise process \mjseqn{\lbrace(i;X_i)\rbrace_{i=0}^{n-1}}
#' with \mjseqn{X_i \sim \mathrm{Normal}(0,\epsilon^2)} is
#' realized and linearly interpolated, yielding a randomly generated
#' function \mjseqn{a(s)} defined continuously on the interval
#' \mjseqn{\lbrack 0,n-1 \rbrack}. A seasonal forcing function
#' `f()` is then defined according to
#'
#' \mjsdeqn{f(s) = \langle\beta\rangle \left\lbrack 1 + \alpha \cos\left(\frac{2 \pi s \Delta t}{\text{365 days}} + a(s)\unicode{x1D7D9}_{\lbrack 0,n-1 \rbrack}\right) \right\rbrack \Delta t}
#'
#' and returned as output. `f()` takes as arguments:
#'
#' \describe{
#'   \item{`s`}{A numeric vector listing values of \mjseqn{s}
#'     at which to evaluate \mjseqn{f(s)}.
#'   }
#'   \item{`par_list`}{A list with numeric scalar elements
#'     `dt_days`, `beta_mean`, and `alpha` listing  values for
#'     \mjseqn{\Delta t}, \mjseqn{\langle\beta\rangle}, and
#'     \mjseqn{\alpha}, respectively.
#'   }
#' }
#'
#' `par_list` does not need an element `epsilon` (giving a value
#' for \mjseqn{\epsilon}), because random number generation and
#' linear interpolation are performed in the enclosing environment
#' of `f()` (the execution environment of `make_beta()`) using the
#' value of `epsilon` found there. That is, `f()` does not need to
#' know `epsilon`, because it does no sampling or interpolating of
#' its own and instead uses the definition of the interpolant
#' \mjseqn{a(s)} that it finds in its enclosing environment.
#'
#' The upshot is that the output of `make_beta()` (i.e., the function
#' `f()`) is randomly generated and reproducible with [base::set.seed()],
#' while the output of `f()` is determined and reproducible without
#' [base::set.seed()]).
#'
#' @param epsilon A numeric scalar. The standard deviation of
#'   the noise process (see Details).
#'
#' @param n An integer scalar. The number of observations in
#'   the noise process (see Details).
#'
#' @return
#' A function with arguments `s` and `par_list` (see Details).
#'
#' @examples
#' epsilon <- pi / 4
#' n <- 1e03
#' set.seed(1734)
#' f <- make_beta(epsilon, n)
#' formals(f)
#' body(f)
#' names(as.list(environment(f)))
#' a <- get("a", envir = environment(f))
#' s <- 0:(n - 1)
#' plot(s, a(s), type = "l", col = "blue")
#'
#' @seealso [make_data()]
#' @export
#' @importFrom stats approxfun rnorm
make_beta <- function(epsilon, n) {
  n <- floor(n)
  a <- approxfun(x = 0:(n - 1), y = rnorm(n, mean = 0, sd = epsilon))
  function(s, par_list) {
    a_val <- ifelse(s < 0 | s > n - 1, 0, a(s))
    one_year <- 365 / par_list$dt_days
    par_list$beta_mean *
      (1 + par_list$alpha * cos(2 * pi * s / one_year + a_val))
  }
}
