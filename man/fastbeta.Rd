% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastbeta.R
\docType{package}
\name{fastbeta}
\alias{fastbeta}
\title{The \pkg{fastbeta} package}
\description{
The \pkg{fastbeta} package
}
\details{
\pkg{fastbeta} implements the FC, S, and SI methods for estimating
time-varying transmission rates
\ifelse{latex}{\out{$\beta(t)$}}{\ifelse{html}{\out{<i>&beta;</i>(<i>t</i>)}}{beta(t)}}
from time series data:
\itemize{
\item \code{\link[=estimate_beta_FC]{estimate_beta_FC()}}
\item \code{\link[=estimate_beta_S]{estimate_beta_S()}}
\item \code{\link[=estimate_beta_SI]{estimate_beta_SI()}}
}

The SI method is substantially more accurate and robust than the
FC and S methods, and should be preferred in practice.

\pkg{fastbeta} additionally implements peak-to-peak iteration (PTPI),
a method for estimating the initial number of susceptibles
\ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}}
from time series data:
\itemize{
\item \code{\link[=ptpi]{ptpi()}}
\item \code{\link[=get_peak_times]{get_peak_times()}}
}

PTPI can be used in conjunction with the FC, S, and SI methods,
which require users to specify an estimate of
\ifelse{latex}{\out{$S_0$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub>}}{S_0}}.

\pkg{fastbeta} includes functions useful for simulating incidence
time series with an underlying, seasonally forced transmission rate.
These may be useful in testing:
\itemize{
\item \code{\link[=make_par_list]{make_par_list()}}
\item \code{\link[=make_data]{make_data()}}
\item \code{\link[=compute_rrmse]{compute_rrmse()}}
}

All methods are based on the SIR model with time-varying rates of
birth, death, and transmission:

\ifelse{latex}{
  \out{
    \begin{array}[rlc]
      S' & = & \nu(t) \widehat{N}_0 - \beta(t) S I - \mu(t) S \\
      I' & = & \beta(t) S I - \gamma I - \mu(t) I \\
      R' & = & \gamma I - \mu(t) R \\
    \end{array}
  }
}{
  \ifelse{html}{
    \out{
      <i>S</i>&prime; = <i>&nu;</i>(<i>t</i>)<i>&Ntilde;</i><sub>0</sub> &minus; <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&mu;</i>(<i>t</i>)<i>S</i><br>
      <i>I</i>&prime; = <i>&beta;</i>(<i>t</i>)<i>SI</i> &minus; <i>&gamma;I</i> &minus; <i>&mu;</i>(<i>t</i>)<i>I</i><br>
      <i>R</i>&prime; = <i>&gamma;I</i> &minus; <i>&mu;</i>(<i>t</i>)<i>R</i>
    }
  }{
    % S' = nu(t)*hatN_0 - beta(t)*S*I - mu(t)*S \cr
    % I' = beta(t)*S*I - gamma*I - mu(t)*I \cr
    % R' = gamma*I - mu(t)*R
  }
}

They are fully described in the referenced manuscript.
}
\references{
deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
time-varying transmission rates for infectious diseases.
}
