#' \loadmathjax
#' The \pkg{fastbeta} package
#'
#' @description
#' \pkg{fastbeta} implements the FC, S, and SI methods for estimating
#' time-varying infectious disease transmission rates \mjseqn{\beta(t)}
#' from disease incidence and mortality data:
#'
#' * [estimate_beta_fc()]
#' * [estimate_beta_s()]
#' * [estimate_beta_si()]
#'
#' The SI method is substantially more robust than the FC and S methods
#' and should be preferred in practice.
#'
#' \pkg{fastbeta} additionally implements peak-to-peak iteration
#' (PTPI), a method for estimating the initial number of susceptibles
#' \mjseqn{S_0} from time series data:
#'
#' * [ptpi()]
#' * [get_peak_times()]
#'
#' PTPI can be used in conjunction with the SI method, which requires
#' users to specify an estimate of \mjseqn{S_0}.
#'
#' \pkg{fastbeta} includes functions useful for simulating incidence
#' time series with an underlying, seasonally forced transmission rate.
#' These may be helpful in testing:
#'
#' * [make_par_list()]
#' * [make_data()]
#' * [compute_rrmse()]
#'
#' All methods are based on the SIR model with time-varying rates of
#' birth, death, and transmission:
#'
#' \mjsdeqn{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu(t) \widehat{N}_0 - \beta(t) S I - \mu(t) S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu(t) I \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R \end{align}}
#'
#' See References for details.
#'
#' @references
#' \insertRef{deJo+20}{fastbeta}
#'
#' @docType package
#' @name fastbeta
#' @importFrom Rdpack reprompt
#' @importFrom mathjaxr preview_rd
NULL
