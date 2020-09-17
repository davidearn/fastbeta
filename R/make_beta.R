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
#' realized and linearly interpolated, yielding a randomly
#' generated function \mjseqn{a(s)} defined continuously on
#' the interval \mjseqn{\lbrack 0,n-1 \rbrack}. A seasonal
#' forcing function `f()` is then defined according to
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
#'     `dt_days`, `beta_mean`, and `alpha` listing values for
#'     \mjseqn{\Delta t}, \mjseqn{\langle\beta\rangle \Delta t},
#'     and \mjseqn{\alpha}, respectively.
#'   }
#' }
#'
#' `par_list` does not need an element `epsilon` (defining \mjseqn{\epsilon}),
#' because `f()` uses the definition of function `a()` (defining \mjseqn{a(s)})
#' that it finds in `environment(f)`. The upshot is that the output of
#' `make_beta()` (i.e., the function `f()`) is randomly generated and
#' reproducible with [base::set.seed()], while the output of `f()` is
#' determined and reproducible without [base::set.seed()].
#'
#' @param epsilon A numeric scalar. Standard deviation of the noise process.
#' @param n An integer scalar. Number of observations in the noise process.
#'
#' @return
#' A closure with arguments `s` and `par_list` (see Details).
#'
#' @examples
#' epsilon <- 0.8
#' n <- 1000
#' pl <- list(dt_days = 7, beta_mean = 1e-05, alpha = 0.05)
#' f <- make_beta(epsilon, n)
#' names(as.list(environment(f)))
#' a <- get("a", envir = environment(f))
#' s <- 0:(n-1)
#' op <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
#' plot(s, a(s), type = "l")
#' plot(s, f(s, pl), type = "l")
#' par(op)
#'
#' @seealso [make_data()]
#' @export
#' @importFrom stats approxfun rnorm
make_beta <- function(epsilon, n) {
  n <- floor(n)
  a <- approxfun(x = 0:(n-1), y = rnorm(n, mean = 0, sd = epsilon))
  function(s, par_list) {
    a_val <- ifelse(s < 0 | s > n-1, 0, a(s))
    one_year <- 365 / par_list$dt_days
    par_list$beta_mean *
      (1 + par_list$alpha * cos(2 * pi * s / one_year + a_val))
  }
}
