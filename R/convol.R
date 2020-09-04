#' Convolve an incidence time series with a delay distribution
#'
#' @description
#' Convolves an equally spaced incidence time series with a discrete
#' distribution of the number of days from infection to reporting,
#' generating an expected reported incidence time series. In addition,
#' simulates reported incidence by sampling from the delay distribution.
#'
#' @details
#' The algorithm does not explicitly require daily observations. For
#' example, `x` can be described most generally as listing the number
#' of infections during each observation interval from interval 0 to
#' interval `length(x)-1`. A different time step (e.g., one week instead
#' of one day) can be chosen as long as the choice is consistent between
#' `x` and `delay_dist`. The output must be interpreted accordingly.
#'
#' @param x An integer vector defining an equally spaced incidence time series.
#'   Lists the number of infections on each day from day 0 to day `length(x)-1`.
#' @param delay_dist A numeric vector. The distribution of the number of days
#'   from infection to reporting. `delay_dist[i]` is the probability that
#'   an infection is reported after `i-1` days. `delay_dist` is replaced with
#'   `delay_dist / sum(delay_dist)` in the event that `sum(delay_dist) != 1`.
#' @param n A non-negative integer scalar. The number of simulations desired.
#'
#' @return
#' A convol object. A list with elements:
#'
#' \describe{
#'   \item{`times`}{A numeric vector equal to `0:(length(x)-1+b)`,
#'     where `b = max(which(delay_dist > 0)) - 1` is the maximum number
#'     of days from infection to reporting, giving time points in days
#'     for the convolution of `x` and `delay_dist`.
#'   }
#'   \item{`x_pad`}{A numeric vector equal to `c(x, rep(NA, b))`.
#'     A convenience allowing one to plot the incidence time series
#'     with `lines(times, x_pad, ...)`.
#'   }
#'   \item{`convolution`}{A numeric vector of length `length(times)`
#'     giving the convolution of `x` and `delay_dist`. `convolution[i]`
#'     is the expected number of infections reported on each day `times[i]`.
#'   }
#'   \item{`simulation`}{If `n = 0`, then `NULL`. Otherwise,
#'     a numeric matrix with `length(times)` rows and `n` columns.
#'     Each column is the result of sampling `delay_dist` `sum(x)` times
#'     and binning the `sum(x)` infections counted in `x` forward in time
#'     according to the sampled delays. `simulation[i, j]` is the number
#'     of infections reported on day `times[i]` in simulation `j` of `n`.
#'   }
#'   \item{`call`}{The function call. The convol object is reproducible
#'     with `eval(call)`.
#'   }
#'   \item{`arg_list`}{A list of the arguments in the function call. The
#'     convol object is reproducible with `do.call(convol, arg_list)`.
#'   }
#' }
#'
#' @examples
#' x <- floor(400 * exp(-seq(-2, 2, by = 0.1)^2))
#' delay_dist <- rep(1 / 11, 11)
#' n <- 5
#' convol_out <- convol(x, delay_dist, n)
#' plot(convol_out)
#'
#' @seealso [methods for class "convol"][convol-methods], [deconvol()]
#' @export
convol <- function(x, delay_dist = c(1), n = 0L) {
  ## Save arguments in a list
  arg_list <- as.list(environment())

  ## Normalize `delay_dist`
  delay_dist <- delay_dist / sum(delay_dist)

  ## Infections occur at times `0:(N-1)`
  N <- length(x)

  ## Time from infection to reporting is at most `b` days
  b <- max(which(delay_dist > 0)) - 1

  ## Pad `x` and `delay_dist` with zeros to length `N + b`
  ## (convenient for matrix operations)
  x <- c(x, rep(0, b))
  delay_dist <- c(delay_dist[1:(b+1)], rep(0, N - 1))

  ## Create a lower triangular transition matrix satisfying
  ## `L[i,j] == delay_dist[i-j+1]`
  L <- array(rep(c(delay_dist, 0), N + b), dim = rep(N + b, 2))
  L[upper.tri(L, diag = FALSE)] <- 0

  ## Convolve `x` and `delay_dist`
  expected_rep_inc <- L %*% x

  if (n > 0) {
    ## Matrix with row `i` listing reporting delays for infection `i`
    delays <- replicate(n,
      sample(
        x       = 0:b,
        size    = sum(x),
        replace = TRUE,
        prob    = delay_dist[1:(b+1)]
      )
    )
    ## Matrix with row `i` listing reporting times for infection `i`
    ## as integers between 1 and `N + b`
    report_times <- rep(seq_along(x), times = x) + delays

    ## Matrix with row `i` listing counts of reports on day `i-1`
    simulated_rep_inc <- apply(report_times, 2,
      function(x) {
        y <- rep(0, N + b)
        y[sort(unique(x))] <- table(x)
        y
      }
    )
  }

  out <- list(
    times       = 0:(N-1+b),
    x_pad       = c(x[1:N], rep(NA, b)),
    convolution = expected_rep_inc,
    simulation  = if (n > 0) simulated_rep_inc else NULL,
    call        = match.call(),
    arg_list    = arg_list
  )
  structure(out, class = c("convol", "list"))
}

