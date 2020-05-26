#' \loadmathjax
#' The \pkg{fastbeta} package
#'
#' @details
#' \pkg{fastbeta} implements the FC, S, and SI methods for estimating
#' time-varying transmission rates \mjeqn{\beta(t)}{beta(t)} from time
#' series data:
#'
#' * [estimate_beta_FC()]
#' * [estimate_beta_S()]
#' * [estimate_beta_SI()]
#'
#' The SI method is substantially more robust than the FC and S methods
#' and should be preferred in practice.
#'
#' \pkg{fastbeta} additionally implements peak-to-peak iteration (PTPI),
#' a method for estimating the initial number of susceptibles
#' \mjeqn{S_0}{S_0} from time series data:
#'
#' * [ptpi()]
#' * [get_peak_times()]
#'
#' PTPI can be used in conjunction with the FC, S, and SI methods,
#' which require users to specify an estimate of \mjeqn{S_0}{S_0}.
#'
#' \pkg{fastbeta} includes functions useful for simulating incidence
#' time series with an underlying, seasonally forced transmission rate.
#' These may be useful in testing:
#'
#' * [make_par_list()]
#' * [make_data()]
#' * [compute_rrmse()]
#'
#' All methods are based on the SIR model with time-varying rates of
#' birth, death, and transmission:
#'
#' \mjdeqn{\begin{align} \frac{\text{d}S}{\text{d}t} &= \nu(t) - \beta(t) S I - \mu(t) S \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu(t) I \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) \end{align}}{(1) dS/dt = nu(t) - beta(t) S I - mu(t) S    (2) dI/dt = beta(t) S I - gamma I - mu(t) I    (3) dR/dt = gamma I - mu(t) R}
#'
#' They are fully described in the referenced manuscript.
#'
#' @references
#' deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
#' time-varying transmission rates for infectious diseases.
#'
#' @md
#' @docType package
#' @name fastbeta
NULL
