test_s2inpars <- function(pars_to_vary   = c("tgen"),
                          par_list_ref   = make_par_list(),
                          scale_factors  = c(1),
                          with_dem_stoch = FALSE,
                          nsim           = 4L,
                          loess_par      = c(NA, NA),
                          ptpi_iter      = 0L) {

## Set-up ==============================================================

# Return to main working directory when you are done
main_wd <- getwd()
on.exit(setwd(main_wd))

# Create a directory for `.RData`, named by simulation method
dirname <- paste0(
  main_wd, "/RData/s2inpars/",
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

# Preallocate memory for output
rrmse <- array(NA,
  dim      = c(dim(par_vals_all), nsim, 2),
  dimnames = list(NULL, colnames(par_vals_all), NULL, c("S", "SI"))
)


## Sensitivity analysis ================================================

for (par in colnames(par_vals_all)) { # loop over parameters

  for (i in 1:nrow(par_vals_all)) { # loop over parameter values

    ## Set-up ----------------------------------------------------------

    # List of values for data-generating parameters.
    # Use reference values.
    par_list_dg <- par_list_ref

    # List of values for input parameters.
    # Use reference values, except for `par`.
    par_list_in <- par_list_ref
    par_list_in[[par]] <- par_vals_all[i, par]

    ## NOTE: This introduces input error, because the input
    ##       value `par_vals_all[i, par]` will differ from
    ##       the data-generating value `par_list_ref[[par]]`
    ##       by a factor of `scale_factors[i]`.

    # Logical matrix: Have you already warned about
    # `estimate_beta_S()` or `estimate_beta_SI()` returning
    # a negative estimate of susceptibles or infecteds?
    been_neg <- array(FALSE,
      dim      = c(2, 2),
      dimnames = list(c("S", "I"), c("S", "SI"))
    )


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

      ## 2. Update value of `S0` in input ... --------------------------
      ##    using peak-to-peak iteration

      if (par == "S0" && ptpi_iter > 0) {
        Z <- estimate_beta_SI(df, par_list_in)$Z
        peaks <- get_peak_times(
          x         = Z,
          period    = with(par_list_in, (365 / 7) / dt_weeks),
          bw_mavg   = 6,
          bw_peakid = 8
        )
        ptpi_out <- ptpi(
          df             = data.frame(Z),
          par_list       = par_list_in,
          a              = with(peaks, phase[1]),
          b              = with(peaks, phase[length(phase)]),
          initial_S0_est = par_list_in$S0,
          iter           = ptpi_iter
        )
        par_list_in$S0 <- ptpi_out$S0_final
      }

      ## 3. Estimate transmission rate ---------------------------------

      estimate_beta <- list(
        S  = estimate_beta_S,
        SI = estimate_beta_SI
      )
      df_est <- suppressWarnings(
        lapply(estimate_beta, function(f) f(df, par_list_in))
      )

      # Logical matrix: Did `estimate_beta_S()` or `estimate_beta_SI()`
      # return a negative estimate of susceptibles or infecteds?
      is_neg <- sapply(df_est,
        function(x) {
          sapply(x[c("S", "I")],
            function(y) {
              any(y < 0, na.rm = TRUE)
            }
          )
        }
      )

      # Warn if an estimate of susceptibles or infecteds is negative
      # and this hasn't already happened for an earlier `k`
      apply(
        expand.grid(c("S", "I"), c("S", "SI"), stringsAsFactors = FALSE), 1,
        function(x) {
          if (is_neg[x[1], x[2]] && !been_neg[x[1], x[2]]) {
            warning(
              par, "_log2f-",
              sprintf("%+05.0f", log(scale_factors[i], 2) * 1000),
              " ... ",
              x[2], " method returns `", x[1], "[i]` < 0 for some `i`",
              call. = FALSE
            )
          }
        }
      )

      # Update `been_neg` to reflect `is_neg`
      been_neg <- been_neg | is_neg

      ## 4. Fit loess curve to raw estimate ----------------------------

      # `loess()` handles `NA` but complains about `NaN` and `Inf`
      df_est <- lapply(df_est,
        function(x) {
          x$beta[is.nan(x$beta) | is.infinite(x$beta)] <- NA
          x
        }
      )

      # An unfortunate quirk in the implementation of `make_data()`
      # and `estimate_beta()`, which could be fixed with much effort:
      # When `trep` is underestimated, `estimate_beta()` fails and the
      # estimate `beta` consists entirely of `NA`. `loess()` throws an
      # error when it encounters this input, so we resort to `try()`.
      loess_fit <- mapply(
        function(x, y) {
          try(
            {
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
            silent = TRUE
          )
        },
        df_est, loess_par,
        SIMPLIFY = FALSE
      )

      ## 4. Compute RRMSE in loess estimate ----------------------------

      # Retain `NA` if susceptibles or infecteds were estimated
      # to be negative or if `loess()` failed earlier
      rrmse[i, par, k, ] <- mapply(
        function(x, y) {
          if (x || inherits(y, "try-error")) {
            return(NA)
          }
          compute_rrmse(df$beta, stats::predict(y))
        },
        apply(is_neg, 2, any), loess_fit
      )

    }

  }

}


## Output ==============================================================

attr(rrmse, "arg_list") <-
  as.list(environment())[names(formals(test_s2inpars))]
rrmse
}
