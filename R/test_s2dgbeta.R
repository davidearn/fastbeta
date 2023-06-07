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
