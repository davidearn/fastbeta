\name{ptpi}
\alias{ptpi}
\title{Estimate
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
using PTPI}
\usage{
ptpi(df = data.frame(), par_list = list(), a, b, initial_S0_est, iter = 0L)
}
\arguments{
\item{df}{A data frame with numeric columns:

\describe{
\item{\code{Z}}{Incidence. \code{Z[i]} is the number of infections between
time
\ifelse{latex}{\out{$t = t_{i-2}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;2</sub>}}{t = t_i-2}}
and time
\ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}.
\code{Z} must be roughly periodic.
}
\item{\code{B}}{Births. \code{B[i]} is the number of births between time
\ifelse{latex}{\out{$t = t_{i-2}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;2</sub>}}{t = t_i-2}}
and time
\ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}.
}
\item{\code{mu}}{Natural mortality rate. \code{mu[i]} is the rate at time
\ifelse{latex}{\out{$t = t_{i-1}$}}{\ifelse{html}{\out{<i>t</i> = <i>t</i><sub><i>i</i>&minus;1</sub>}}{t = t_i-1}}
expressed per unit
\ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
and per capita.
}
}

\code{B} is optional if \code{hatN0} and \code{nu} are defined in \code{par_list}, and
\code{mu} is optional if \code{mu} is defined in \code{par_list} (see Details).}

\item{par_list}{A list containing:

\describe{
\item{\code{hatN0}}{[ \ifelse{latex}{\out{$\widehat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}} ]
Population size at time
\ifelse{latex}{\out{$t = 0$}}{\ifelse{html}{\out{<i>t</i> = 0}}{t = 0}}
years.
}
\item{\code{nu}}{[ \ifelse{latex}{\out{$\nu_\text{c}$}}{\ifelse{html}{\out{<i>&nu;<sub>c</sub></i>}}{nu_c}} ]
Birth rate expressed per unit
\ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
and relative to
\ifelse{latex}{\out{$\hat{N}_0$}}{\ifelse{html}{\out{<i>&Ntilde;</i><sub>0</sub>}}{hatN_0}}
(if modeled as constant).
}
\item{\code{mu}}{[ \ifelse{latex}{\out{$\mu_\text{c}$}}{\ifelse{html}{\out{<i>&mu;</i><sub>c</sub>}}{mu_c}} ]
Natural mortality rate expressed per unit
\ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
and per capita (if modeled as constant).
}
}

\code{hatN0} and \code{nu} are optional if \code{B} is defined in \code{df}, and
\code{mu} is optional if \code{mu} is defined in \code{df} (see Details).}

\item{a}{Integer scalar. Index of first peak in \code{df$Z}, possibly
obtained via \code{get_peak_times()}.}

\item{b}{Integer scalar. Index of last peak in \code{df$Z} in phase
with first peak, possibly obtained via \code{get_peak_times()}.}

\item{initial_S0_est}{Numeric scalar. An initial estimate of
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.}

\item{iter}{Integer scalar. The number of estimates of
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
to generate before stopping.}
}
\value{
A list containing:

\describe{
\item{\code{S_mat}}{A numeric matrix with dimensions
\code{c(nrow(df), iter+1)} containing the susceptible time series
generated in each iteration.
}
\item{\code{S0}}{A numeric vector listing in order all \code{1+iter}
estimates of
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.
Equal to \code{S_mat[1, ]}.
}
\item{\code{S0_final}}{The final estimate of
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}.
Equal to \code{S0[length(S0)]}.
}
\item{\code{SA}}{A numeric vector listing in order all \code{1+iter}
estimates of
\ifelse{latex}{\out{$S_\text{a} = S(t_\text{a})$}}{\ifelse{html}{\out{<i>S</i><sub>a</sub> = <i>S</i>(<i>t</i><sub>a</sub>)}}{S_a = S(t_a)}},
where
\ifelse{latex}{\out{$t_\text{a}$}}{\ifelse{html}{\out{<i>t</i><sub>a</sub>}}{t_a}}
is the time point corresponding to row \code{a} in \code{df}.
Equal to \code{S_mat[a, ]}.
}
\item{\code{SA_final}}{The final estimate of
\ifelse{latex}{\out{$S_\text{a} = S(t_\text{a})$}}{\ifelse{html}{\out{<i>S</i><sub>a</sub> = <i>S</i>(<i>t</i><sub>a</sub>)}}{S_a = S(t_a)}}.
Equal to \code{SA[length(SA)]}.
}
}

A list of the arguments of \code{ptpi()} is included as an attribute.
}
\description{
Using the method of peak-to-peak iteration (PTPI, see References),
\code{ptpi()} estimates the initial number of susceptibles
\ifelse{latex}{\out{$S_0 = S(t_0)$}}{\ifelse{html}{\out{<i>S</i><sub>0</sub> = <i>S</i>(<i>t</i><sub>0</sub>)}}{S_0 = S(t_0)}}
from time series of incidence (\emph{must be roughly periodic}), births,
and natural mortality, observed at equally spaced time points
\ifelse{latex}{\out{$t_k = t_0 + k \Delta t$}}{\ifelse{html}{\out{<i>t<sub>k</sub></i> = <i>t</i><sub>0</sub>+<i>k&Delta;t</i>}}{t_k = t_0 + k*Dt}}
(for \ifelse{latex}{\out{$k = 0,\ldots,n$}}{\ifelse{html}{\out{<i>k</i> = 0,...,<i>n</i>}}{k = 0,...,n}}),
where
\ifelse{latex}{\out{$\Delta t$}}{\ifelse{html}{\out{<i>&Delta;t</i>}}{Dt}}
denotes the observation interval.
}
\section{Mock vital data}{

If \code{df$B} is undefined in the function call, then \code{df$B[i]} gets the
value \code{with(par_list, nu * hatN0 * 1)} for all \code{i}. If \code{df$mu} is
undefined in the function call, then \code{df$mu[i]} gets the value
\code{with(par_list, mu)} for all \code{i}.
}

\section{Missing data}{

Missing values in \code{df[, c("Z", "B", "mu")]} are mostly not tolerated.
At the moment, \code{ptpi()} makes no effort to impute them, so imputation
must be done beforehand.
}

\examples{
# Simulate 20 years of disease incidence,
# observed weekly
par_list <- make_par_list(dt_weeks = 1)
df <- make_data(
  par_list = par_list,
  n = 20 * 365 / 7, # 20 years is ~1042 weeks
  with_dem_stoch = TRUE,
  seed = 5
)

# Plot incidence time series, and note the
# apparent 1-year period
plot(Z ~ t_years, df,
  type = "l",
  xlab = "Time (years)",
  ylab = "Incidence"
)

# Find peaks in incidence time series
peaks <- get_peak_times(
  x = df$Z,
  period = 365 / 7, # 1 year is ~52 weeks
  bw_mavg = 6,
  bw_peakid = 8
)

# Verify that peaks were identified correctly
abline(v = df$t_years[peaks$all], lty = 2, col = "red")

# Index of first peak
a <- with(peaks, phase[1])

# Index of last peak in phase with first
b <- with(peaks, phase[length(phase)])

# Estimate `S0` from incidence via PTPI,
# starting from an erroneous initial guess
# (mock vital data generated internally)
ptpi_out <- ptpi(
  df = df["Z"],
  par_list = par_list,
  a = a,
  b = b,
  initial_S0_est <- df$S[1] * 4,
  iter = 25
)

# Sequence of estimates
ptpi_out$S0

# Relative error in final estimate
(ptpi_out$S0_final - df$S[1]) / df$S[1]

}
\references{
deJonge MS, Jagan M, Krylova O, Earn DJD. Fast estimation of
time-varying transmission rates for infectious diseases.
}
