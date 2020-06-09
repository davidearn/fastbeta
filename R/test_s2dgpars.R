#' Measure sensitivity to data-generating parameters
#'
#' `test_s2dgpars()` measures the sensitivity of transmission rate
#' estimation error to parameters used to simulate reported incidence
#' data. It varies the data-generating parameters *one at a time*, and,
#' using each parametrization, it
#' (i) performs simulations,
#' (ii) estimates the data-generating, seasonally forced transmission
#' rate from the simulated data, without input error, and
#' (iii) records the estimation error.
#' Simulations are saved as `.RData` in the directory
#' `"./RData/s2dgpars/"`.
#'
#' @details
#' The space of data-generating parameters is sampled by assigning a
#' parameter in `pars_to_vary` its reference value in `par_list_ref`
#' multiplied by a scale factor in `scale_factors`, and assigning all
#' other parameters their reference value in `par_list_ref`. The total
#' number of parametrizations considered is thus given by
#' `length(pars_to_vary) * length(scale_factors)`. For each
#' parametrization
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
#' There is one caveat in the above description: as `nu`, `mu`, and
#' `tgen` vary, values for `Rnaught`, `N0`, `S0` and `I0` are *not*
#' not taken from `par_list_ref`. As
#' \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}},
#' \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}},
#' and
#' \ifelse{latex}{\out{$t_\text{gen}$}}{\ifelse{html}{\out{<i>t</i><sub>gen</sub>}}{t_gen}}
#' vary,
#' \ifelse{latex}{\out{$\mathcal{R}_0$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub>}}{calR_0}}
#' varies concurrently in order to enforce the identity
#'
#' \ifelse{latex}{\out{$\mathcal{R}_0 = \frac{\nu_\text{c} \widehat{N}_0}{\mu_\text{c}} \cdot \frac{\langle\beta\rangle}{\gamma + \mu}\,.$}}{\ifelse{html}{\out{<i>&Rscr;</i><sub>0</sub> = (<i>&nu;</i><sub>c</sub> <i>&Ntilde;</i> / <i>&mu;</i><sub>c</sub>)(&langle;<i>&beta;</i>&rangle; / (<i>&gamma;</i> + <i>&mu;</i><sub>c</sub>)).}}{calR_0 = ((nu_c*hatN0)/mu_c)*(<beta>/(gamma + mu)).}}
#'
#' Similarly,
#' \ifelse{latex}{\out{$N_0$}}{\ifelse{html}{\out{<i>N</i><sub>0</sub>}}{N_0}},
#' \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}},
#' and
#' \ifelse{latex}{\out{$I_0$}}{\ifelse{html}{\out{<i>I</i><sub>0</sub>}}{I_0}}
#' vary concurrently in order to ensure that they specify a point very
#' near the attractor of the system of SIR equations. (The system and its
#' attractor change as functions of
#' \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}},
#' \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}},
#' and
#' \ifelse{latex}{\out{$t_\text{gen}$}}{\ifelse{html}{\out{<i>t</i><sub>gen</sub>}}{t_gen}}.)
#'
#' @param pars_to_vary Character vector. Contains the names of all
#'   data-generating parameters to be varied, chosen from `"S0"`,
#'   `"I0"`, `"nu"`, `"mu"`, `"tgen"`, and `"prep"`.
#' @param par_list_ref A list returned by
#'   `make_par_list(..., beta_mean = NA, N0 = NA, S0 = NA, I0 = NA)`.
#'   Contains reference values for all data-generating parameters.
#'   The values listed for `Rnaught`, `N0`, `S0`, and `I0` are not
#'   used when `nu`, `mu`, or `tgen` acts as the independent variable
#'   in the analysis (see Details).
#' @param scale_factors Numeric vector. The values of `pars_to_vary[j]`
#'   considered in the analysis are
#'   `par_list_ref[[pars_to_vary[j]]] * scale_factors`.
#' @param with_dem_stoch Logical scalar, passed to [make_data()]. If
#'   `TRUE`, then simulations account for demographic stochasticity.
#' @param nsim Integer scalar. The number of simulations to perform
#'   using each parametrization.
#' @param loess_par Integer vector of length 2. Determines the degree
#'   of smoothing when loess curves are fit to raw estimates of the
#'   seasonally forced transmission rate (see Details). Estimates
#'   generated by [estimate_beta_S()] use `loess_par[1]`, while those
#'   generated by [estimate_beta_SI()] use `loess_par[2]`.
#' @param only_make_data Logical scalar. If `TRUE`, then simulations
#'   are performed and saved, as usual, but nothing more is done and
#'   nothing is returned.
#'
#' @return
#' Unless `only_make_data` is `TRUE`, a numeric array with dimensions
#' `c(length(scale_factors), length(pars_to_vary), nsim, 2)`. Stores
#' the relative root mean square error (RRMSE) in each estimate of the
#' data-generating, seasonally forced transmission rate. The
#' `[i, j, k, m]`th entry corresponds to simulation `k` of `nsim` using
#' `par_vals_all[i, j]`, and estimate `m` of 2 from that simulation
#' (S method for `m = 1`, SI method for `m = 2`).
#'
#' A list of the arguments of `test_s2dgpars()` is included as an
#' attribute.
#'
#' @md
#' @export
test_s2dgpars <- function(pars_to_vary   = c("tgen"),
                          par_list_ref   = make_par_list(),
                          scale_factors  = c(1),
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
  main_wd, "/RData/s2dgpars/",
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

# Generate the set of values considered for each parameter
# listed in `pars_to_vary`
par_vals_ref <- unlist(par_list_ref[pars_to_vary])
par_vals_all <- as.data.frame(outer(scale_factors, par_vals_ref))

# Preallocate memory for output, unless you only want simulations
if (!only_make_data) {
  rrmse <- array(NA,
    dim      = c(dim(par_vals_all), nsim, 2),
    dimnames = list(NULL, colnames(par_vals_all), NULL, c("S", "SI"))
  )
}


## Sensitivity analysis ================================================

for (par in colnames(par_vals_all)) { # loop over parameters

  for (i in 1:nrow(par_vals_all)) { # loop over parameter values

    ## Set-up ----------------------------------------------------------

    # Create a subdirectory for `.RData`, named by parametrization
    subdirname <- paste0(
      dirname,
      # current parameter
      par,
      # log2 ratio of current value to reference value
      "_log2f-", sprintf("%+05.0f", log(scale_factors[i], 2) * 1000),
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

      # Update `par`
      mpl_args[[par]] <- par_vals_all[i, par]

      # If you just updated `nu`, `mu`, or `tgen`, then the system
      # being simulated has changed. This means that:
      # * `Rnaught` in `par_list_ref` is not the correct basic
      #   reproduction number for the new system, assuming that
      #   `beta_mean` is the same as in `par_list_ref`.
      # * `N0`, `S0`, and `I0` in `par_list_ref` do not specify
      #   a point near the attractor of the new system.
      # In this case, erase them from `mpl_args`. This will instruct
      # `make_par_list()` to compute and assign the desired values.
      if (par %in% c("nu", "mu", "tgen")) {
        mpl_args[["Rnaught"]] <- NA
        mpl_args[c("N0", "S0", "I0")] <- NA
      }
      
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
        "par ", which(colnames(par_vals_all) == par), " of ",
        ncol(par_vals_all), ", ",
        "val ", i, " of ", nrow(par_vals_all), ", ",
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

      rrmse[i, par, k, ] <- sapply(loess_fit,
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
  as.list(environment())[names(formals(test_s2dgpars))]
rrmse
}
