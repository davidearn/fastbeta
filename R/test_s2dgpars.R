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
