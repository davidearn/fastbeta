#' The \pkg{fastbeta} package
#'
#' @details
#' \pkg{fastbeta} implements the FC, S, and SI methods for estimating
#' time-varying transmission rates
#' \ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
#' from time series data:
#'
#' * [estimate_beta_FC()]
#' * [estimate_beta_S()]
#' * [estimate_beta_SI()]
#'
#' The SI method is substantially more accurate and robust than the
#' FC and S methods, and should be preferred in practice.
#'
#' \pkg{fastbeta} additionally implements peak-to-peak iteration (PTPI),
#' a method for estimating the initial number of susceptibles
#' \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}}
#' from time series data:
#'
#' * [ptpi()]
#' * [get_peak_times()]
#'
#' PTPI can be used in conjunction with the FC, S, and SI methods,
#' which require users to specify an estimate of
#' \ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}}.
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
#' \ifelse{latex}{
#'   \out{
#'     \begin{array}[rlc]
#'       S' & = & \nu(t) \widehat{N}_0 - \beta(t) S I - \mu(t) S \\
#'       I' & = & \beta(t) S I - \gamma I - \mu(t) I \\
#'       R' & = & \gamma I - \mu(t) R \\
#'     \end{array}
#'   }
#' }{
#'   \ifelse{html}{
#'     \out{
#'       <i>S</i>&prime; = <i>&nu;</i>(<i>t</i>)<i>&Ntilde;</i><sub>0</sub> &minus; <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&mu;</i>(<i>t</i>)<i>S</i><br>
#'       <i>I</i>&prime; = <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&gamma;I</i> &minus; <i>&mu;</i>(<i>t</i>)<i>I</i><br>
#'       <i>R</i>&prime; = <i>&gamma;I</i> &minus; <i>&mu;</i>(<i>t</i>)<i>R</i>
#'     }
#'   }{
#'     % S' = nu(t)*hatN_0 - beta(t)*S*I - mu(t)*S \cr
#'     % I' = beta(t)*S*I - gamma*I - mu(t)*I \cr
#'     % R' = gamma*I - mu(t)*R
#'   }
#' }
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
