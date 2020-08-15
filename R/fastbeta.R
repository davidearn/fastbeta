#' Estimate time-varying transmission rates
#'
#' @description
#' Generates a fastbeta object.
#'
#' @details
#' Details to come.
#'
#' @param df A data frame with time series data.
#' @param par_list A list of parameter values.
#'
#' @return
#' A fastbeta object.
#'
#' @export
fastbeta <- function(df, par_list) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")  
  } else if (!all(c("t", "Z", "B", "mu") %in% names(df))) {
    stop("`df` is missing necessary columns.")
  } 
  if (!is.list(par_list)) {
    stop("`par_list` must be a list.")
  } else if (!all(c("S0", "I0", "tgen") %in% names(par_list))) {
    stop("`par_list` is missing necessary elements.")
  }
  
  # Save arguments in a list
  arg_list <- as.list(environment())
    
  # Missing values are not tolerated. Zeros in incidence
  # *are* tolerated but can have undesired effects. See
  # `?estimate_beta_si`.
  df[c("Z", "B", "mu")] <- mapply(impute_na,
    x = df[c("Z", "B", "mu")],
    zero_as_na = c(TRUE, FALSE, FALSE)
  )
  out <- estimate_beta_si(df, par_list)

  # Negative susceptibles indicates that incidence was
  # overestimated or births were underestimated or both
  if (any(out$S < 0, na.rm = TRUE)) {
    warning(
      "Negative elements in susceptibles (`S`) column. Retry with:",
      "\n* scaled down incidence (`df$Z`), and/or",
      "\n* scaled up births (`df$B`)."
    )
  }

  structure(out,
    class = c("fastbeta", "data.frame"),
    call = match.call(),
    arg_list = arg_list
  )
}
