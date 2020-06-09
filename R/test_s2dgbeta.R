#' Measure sensitivity to the data-generating transmission rate
#'
#' `test_s2dgbeta()` measures the sensitivity of transmission rate
#' estimation error to the basic reproduction number
#' \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}
#' and seasonal amplitude
#' \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}}
#' used to simulate reported incidence data. It considers a grid in the
#' \ifelse{latex}{\out{$(\mathcal{R}_0,\alpha)$}}{\ifelse{html}{\out{(<i>&Rscr;</i><sub>0</sub>,<i>&alpha;</i>)}}{(calR_0,alpha)}}
#' parameter subspace. Using each
#' \ifelse{latex}{\out{$(\mathcal{R}_0,\alpha)$}}{\ifelse{html}{\out{(</i>&Rscr;</i><sub>0</sub>,<i>&alpha;</i>)}}{(calR_0,alpha)}}
#' pair, it
#' (i) performs simulations,
#' (ii) estimates the data-generating, seasonally forced transmission
#' rate from the simulated data, without input error, and
#' (iii) records the estimation error.
#' Simulations are saved as `.RData` in the directory
#' `"./RData/s2dgbeta/"`.
#'
#' @details
#' The space of data-generating parameters is sampled by exploring a
#' grid of
#' \ifelse{latex}{\out{$(\mathcal{R}_0,\alpha)$}}{\ifelse{html}{\out{(<i>&Rscr;</i><sub>0</sub>,<i>&alpha;</i>)}}{(calR_0,alpha)}}
#' pairs created from the vectors `Rnaught` and `alpha`, and assigning
#' all other parameters their reference value in `par_list_ref`.
#' The total number of parametrizations considered is thus given by
#' `length(Rnaught) * length(alpha)`. For each parametrization
#' \ifelse{latex}{\out{$\mathbf{\theta}$}}{\ifelse{html}{\out{<b><i>&theta;</i></b>}}{theta}}:
#'
#' 1. [make_data()] is called to simulate `nsim` reported incidence
#'    time series, with arguments:
#'
#'    * `par_list = par_list_dg`, indicating the data-generating
#'      parameter values
#'      \ifelse{latex}{\out{$\mathbf{\theta}$}}{\ifelse{html}{\out{<b><i>&theta;</i></b>}}{theta}}.
#'    * `n = 1042`, indicating the time series length.
#'    * `with_dem_stoch = with_dem_stoch`, indicating whether the
#'      simulation should account for demographic stochasticity.
#'    * `seed = k`, where `k` is the simulation count out of `nsim`,
#'      making the result of each simulation reproducible.
#'
#' 2. [estimate_beta_S()] and [estimate_beta_SI()] are called to
#'    estimate the seasonally forced transmission rate from each
#'    simulation, with arguments:
#'
#'    * `df = df`, where `df` is the output of [make_data()],
#'      supplying the simulated reported incidence time series.
#'    * `par_list = par_list_dg`, indicating the input parameter values
#'      \ifelse{latex}{\out{$\mathbf{\xi}$}}{\ifelse{html}{\out{<b><i>&xi;</i></b>}}{xi}}.
#'      There is no error in the input. In other words, the input
#'      parameter values are identical to the data-generating
#'      parameter values.
#'
#'    Mock (constant) birth and natural mortality time series are
#'    created internally, with `B[i] = with(par_list, hatN0 * nu * 1)`
#'    and `mu[i] = with(par_list, mu)` for all `i`.
#' 3. [stats::loess()] is called to fit a loess curve to each
#'    raw estimate of the seasonally forced transmission rate,
#'    with arguments
#'
#'    * `formula = beta ~ t` and `data = df_est`, where `df_est` is
#'      the output of [estimate_beta_S()] or [estimate_beta_SI()],
#'      indicating the time series to which the loess curve should
#'      be fit.
#'    * `span = loess_par[i] / nrow(df_est)`, indicating (roughly)
#'      the proportion of points to be weighted in local regression.
#'    * `degree = 2`, indicating that a quadratic polynomial should
#'      be fit locally.
#'    * `na.action = "na.exclude"`, indicating that missing values
#'      should be omitted in local regression but preserved in the
#'      output of [stats::predict()].
#'
#' 4. [compute_rrmse()] is called to compute the error in each
#'    loess estimate, with arguments
#'
#'    * `x = df$beta`, indicating the true value of the transmission
#'      rate at each observation time.
#'    * `y = predict(loess_fit)`, where `loess_fit` is the output
#'      of [stats::loess()], indicating the estimated value of the
#'      transmission rate at each observation time.
#'
#' There is one caveat in the above description: values for `beta_mean`,
#' `N0`, `S0`, and `I0` are *not* taken from `par_list_ref`. As
#' \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}
#' varies,
#' \ifelse{latex}{\out{$\langle\beta\rangle$}}{\ifelse{html}{\out{&langle;<i>&beta;</i>&rangle;}}{<beta>}}
#' varies concurrently in order to enforce the identity
#'
#' \ifelse{latex}{\out{$\mathcal{R}_0 = \frac{\nu_\text{c} \widehat{N}_0}{\mu_\text{c}} \cdot \frac{\langle\beta\rangle}{\gamma + \mu}\,.$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub> = (<i>&nu;</i><sub>c</sub> <i>&Ntilde;</i> / <i>&mu;</i><sub>c</sub>)(&langle;<i>&beta;</i>&rangle; / (<i>&gamma;</i> + <i>&mu;</i><sub>c</sub>)).}}{calR_0 = ((nu_c*hatN0)/mu_c)*(<beta>/(gamma + mu)).}}
#'
#' Similarly, as
#' \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}
#' and
#' \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}}
#' vary,
#' \ifelse{latex}{\out{$N_0$}}{\ifelse{html}{\out{<i>N</i><sub>0</sub>}}{N_0}},
#' \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}},
#' and
#' \ifelse{latex}{\out{$I_0$}}{\ifelse{html}{\out{<i>I</i><sub>0</sub>}}{I_0}}
#' vary in order to ensure that they specify a point very
#' near the attractor of the system of SIR equations.
#' (The system and its attractor change as functions of
#' \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}
#' and
#' \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}}.)
#'
#' @param Rnaught Numeric vector. Values considered for the basic
#'   reproduction number
#'   \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}.
#' @param alpha Numeric vector. Values considered for the seasonal
#'   amplitude
#'   \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<i>&alpha;</i>}}{alpha}}.
#' @param par_list_ref A list returned by `make_par_list(...)`.
#'   Contains values for all data-generating parameters. The values
#'   listed for `Rnaught`, `alpha`, `beta_mean`, `N0`, `S0`, and `I0`
#'   are not used (see Details).
#' @param with_dem_stoch Logical scalar, passed to [make_data()]. If
#'   `TRUE`, then simulations account for demographic stochasticity.
#' @param nsim Integer scalar. The number of simulations to perform
#'   using each parametrization.
#' @param loess_par Integer vector of length 2. Determines the degree
#'   of smoothing when loess curves are fit to raw estimates of the
#'   seasonally forced transmission rate (see Details). Estimates
#'   generated by [estimate_beta_S()] use `loess_par[1]`, while those
#'   generated by [estimate_beta_SI()] use `loess_par[2]`.
#' @param only_make_data Logical. If `TRUE`, then simulations are
#'   performed and saved, as usual, but nothing more is done and
#'   nothing is returned.
#'
#' @return
#' Unless `only_make_data` is `TRUE`, a numeric array with dimensions
#' `c(length(Rnaught), length(alpha), nsim, 2)`. Stores the relative
#' root mean square error (RRMSE) in each estimate of the
#' data-generating, seasonally forced transmission rate.
#' The `[i, j, k, m]`th entry corresponds to simulation `k` of `nsim`
#' using `Rnaught[i]` and `alpha[j]`, and estimate `m` of 2 from that
#' simulation (S method for `m = 1`, SI method for `m = 2`).
#'
#' A list of the arguments of `test_s2dgbeta()` is included as an
#' attribute.
#'
#' @md
#' @export
test_s2dgbeta <- function(Rnaught        = c(20),
                          alpha          = c(0.08),
                          par_list_ref   = make_par_list(),
                          with_dem_stoch = FALSE,
                          nsim           = 4L,
                          loess_par      = c(NA, NA),
                          only_make_data = FALSE) {

## Set-up ==============================================================

# Return to main working directory when you are done
main_wd <- getwd()
on.exit(setwd(main_wd))

# Create a directory for `.RData`, named by simulation method
dirname <- paste0(
  main_wd, "/RData/s2dgbeta/",
  # with or without environmental stochasticity
  if (with(par_list_ref, epsilon == 0)) "xx" else "es",
  # with or without demographic stochasticity
  if (with_dem_stoch) "ds" else "xx",
  # with or without observation error
  if (with(par_list_ref, prep == 1 && trep == 0)) "xx" else "oe",
  "/"
)
if (!dir.exists(dirname)) {
  dir.create(dirname, recursive = TRUE)
}
setwd(dirname)

# Preallocate memory for output, unless you only want simulations
if (!only_make_data) {
  rrmse <- array(NA,
    dim      = c(length(Rnaught), length(alpha), nsim, 2),
    dimnames = list(NULL, NULL, NULL, c("S", "SI"))
  )
}


## Sensitivity analysis ================================================

for (i in seq_along(Rnaught)) {

  for (j in seq_along(alpha)) {

    ## Set-up ----------------------------------------------------------

    # Create a subdirectory for `.RData`, named by parametrization
    subdirname <- paste0(
      dirname,
      "Rnaught-", sprintf("%03.0f", Rnaught[i] * 10),
      "_alpha-", sprintf("%04.0f", alpha[j] * 1000),
      "/"
    )
    if (!dir.exists(subdirname)) {
      dir.create(subdirname, recursive = TRUE)
    }
    setwd(subdirname)

    # Make a list of values for data-generating parameters,
    # if you haven't already. Starting from scratch is costly:
    # the desired initial state is the state after a transient,
    # which takes several seconds to compute.
    if (file.exists("par_list_dg.RData")) {
      load("par_list_dg.RData")
    } else {
      # Construct a list of arguments to pass to `make_par_list()`
      # by modifying `par_list_ref`
      mpl_args <- par_list_ref

      # Update `Rnaught` and `alpha`
      mpl_args["Rnaught"] <- Rnaught[i]
      mpl_args["alpha"] <- alpha[j]

      # Erase `beta_mean`. This instructs `make_par_list()` to
      # compute and assign a value based on the new `Rnaught`.
      mpl_args[["beta_mean"]] <- NA

      # Erase `N0`, `S0`, and `I0`. This instructs `make_par_list()`
      # to compute and assign the state of the new system of SIR
      # equations after a transient.
      mpl_args[c("N0", "S0", "I0")] <- NA

      # Finally, supply the arguments in a call to `make_par_list()`
      par_list_dg <- do.call(make_par_list, mpl_args)

      # Save the output
      save(par_list_dg, file = "par_list_dg.RData")
    }

    # Make a list of values for input parameters. Consider
    # the ideal case in which there is no input error.
    par_list_in <- par_list_dg


    ## Everything else -------------------------------------------------

    for (k in 1:nsim) {

      message(
        "`Rnaught` val ", i, " of ", length(Rnaught), ", ",
        "`alpha` val ", j, " of ", length(alpha), ", ",
        "sim ", k, " of ", nsim
      )

      ## 1. Simulate data ----------------------------------------------

      filename <- paste0("sim", sprintf("%04.0f", k), ".RData")
      if (file.exists(filename)) {
        load(filename)
      } else {
        df <- make_data(
          par_list       = par_list_dg,
          n              = 20 * 365 / 7,
          with_dem_stoch = with_dem_stoch,
          seed           = k
        )
        save(df, file = filename)
      }

      # Stop here if you only want the simulation
      if (only_make_data) {
        next
      }

      ## 2. Estimate transmission rate ---------------------------------

      estimate_beta <- list(
        S  = estimate_beta_S,
        SI = estimate_beta_SI
      )
      df_est <- lapply(estimate_beta, function(f) f(df, par_list_in))

      ## 3. Fit loess curve to raw estimate ----------------------------

      # `loess()` handles `NA` but complains about `NaN` and `Inf`
      df_est <- lapply(df_est,
        function(x) {
          x$beta[is.nan(x$beta) | is.infinite(x$beta)] <- NA
          x
        }
      )

      loess_fit <- mapply(
        function(x, y) {
          stats::loess(
            formula   = beta ~ t,
            data      = x,
            span      = y / nrow(x),
            degree    = 2,
            na.action = "na.exclude",
            control   = stats::loess.control(
              surface    = "direct",
              statistics = "none"
            )
          )
        },
        df_est, loess_par,
        SIMPLIFY = FALSE
      )

      ## 4. Compute RRMSE in loess estimate ----------------------------

      rrmse[i, j, k, ] <- sapply(loess_fit,
        function(x) compute_rrmse(df$beta, stats::predict(x))
      )

    }

  }

}


## Output ==============================================================

if (only_make_data) {
  return(invisible())
}

attr(rrmse, "arg_list") <-
  as.list(environment())[names(formals(test_s2dgbeta))]
rrmse
}
