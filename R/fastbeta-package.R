#' \loadmathjax
#' fastbeta: Fast Estimation of Time-Varying Transmission Rates
#'
#' @description
#' Fast methods for estimating time-varying infectious disease
#' transmission rates from incidence and mortality time series.
#'
#' @details
#' The \pkg{fastbeta} package implements and extends the methods described
#' in \insertCite{Jaga+20;textual}{fastbeta}. The package is modular, with
#' the core machinery divided across different functions, each with its own
#' purpose and convenient plotting methods.
#'
#' [fastbeta()] is a wrapper giving access to the **FC, S, SI, and
#' SEI methods** of estimating time-varying transmission rates from
#' incidence time series. The FC, S, and SI methods are derived from
#' the SIR model
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}I}{\text{d}t} &= \beta(t) S I - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,,\end{align*}}
#'
#' while the SEI method
#' (not discussed in \insertCite{Jaga+20;textual}{fastbeta})
#' is a generalization of the SI method to the SEIR model
#'
#' \mjsdeqn{\begin{align*} \frac{\text{d}S}{\text{d}t} &= \nu(t) - \beta(t) S I - \mu(t) S\,, \cr \frac{\text{d}E}{\text{d}t} &= \beta(t) S I - \sigma E - \mu(t) E\,, \cr \frac{\text{d}I}{\text{d}t} &= \sigma E - \gamma I - \mu(t) I\,, \cr \frac{\text{d}R}{\text{d}t} &= \gamma I - \mu(t) R\,.\end{align*}}
#'
#' Details of the individual algorithms are gathered
#' [here][estimate-beta].
#'
#' As a result of process and observation error, [fastbeta()]
#' can generate noisy transmission rate estimates from which
#' it is difficult to discern temporal patterns of interest.
#' [try_loess()] greatly facilitates the process of **fitting
#' smooth loess curves** to noisy time series using different
#' values for the smoothing parameter.
#'
#' [bootbeta()] can be used to generate **bootstrap confidence intervals**
#' on transmission rate estimates produced by [fastbeta()] (or on loess fits
#' to those estimates).
#'
#' In the typical case where one knows reported incidence but not
#' incidence, [fastbeta()] should not be used directly. Instead,
#' [deconvol()], which implements the **deconvolution** algorithm
#' of \insertCite{Gold+09;textual}{fastbeta}, can be used to
#' reconstruct incidence from reported incidence. The resulting
#' deconvolved incidence time series can then be passed to
#' [fastbeta()]. [convol()] can be used for testing [deconvol()].
#'
#' [ptpi()] implement the **peak-to-peak iteration** introduced in
#' \insertCite{Jaga+20;textual}{fastbeta}. If the initial number
#' of susceptible individuals (a parameter to which [fastbeta()]
#' is sensitive) is not known, and if incidence is roughly periodic,
#' [ptpi()] can be used to produce a reasonable estimate starting
#' from a poor initial guess.
#'
#' [make_data()] can be used to simulate **epidemic time series
#' data** against which many of these functions can be tested.
#' See [fastbeta()], [bootbeta()], and [ptpi()] for examples.
#'
#' @references
#' \insertRef{Jaga+20}{fastbeta}
#'
#' \insertRef{FineClar82}{fastbeta}
#'
#' \insertRef{Gold+09}{fastbeta}
#'
#' @docType package
#' @keywords internal
#' @name fastbeta-package
#' @importFrom Rdpack reprompt
#' @importFrom mathjaxr preview_rd
NULL

## For `R CMD check`
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("tgen", "tlat", "tinf",
      "S", "E", "logE", "logI", "R", "Zcum", "Bcum")
  )
}
