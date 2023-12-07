\name{NEWS}
\title{News for \R{} Package \pkg{fastbeta}}

\section{Changes in version 0.2.0 (2023-12-07)}{
  \itemize{
    \item Three backwards incompatible changes to the API
    due to generalization allowing for loss of immunity:
    \itemize{
      \item \code{fastbeta} and \code{fastbeta.bootstrap}
      expect \code{constants = \code{c(gamma, delta, S0, I0, R0)}},
      no longer \code{constants = \code{c(gamma, S0, I0)}}.

      \item \code{ptpi} gains argument \code{constants} and loses
      argument \code{start}.
      It expects \code{constants = c(gamma, delta, Sa, Ia, Ra)}
      instead of \code{start = Sa}.

      \item \code{fastbeta} returns (modulo attributes)
      \code{cbind(S, I, R, beta)}, no longer \code{cbind(beta, S, I)}.
    }

    \item \file{data/*.R} are no longer copied into \file{inst/scripts}
    at install time by a \command{make} rule in \file{src/Makevars}.
    There is no need: the scripts are preserved due to the setting of
    \code{BuildResaveData: no}.  Hence \file{src/Makevars} is removed.
  }
}

\section{Changes in version 0.1.1 (2023-11-28)}{
  \itemize{
    \item Changes are logged in \file{inst/NEWS.Rd}.

    \item \code{Depends: R (>= 4.3.0)} for \verb{amsmath} in PDF
    output, which was not supported until 4.2.2.

    \item \file{DESCRIPTION} sets \code{BuildResaveData: no},
    so that \command{R CMD build} works without
    \option{--no-resave-data} when \pkg{fastbeta} is not installed.

    \item \file{data/*.R} are copied into \file{inst/scripts}
    at install time by a \command{make} rule in \file{src/Makevars}.

    \item Typo in \file{src/sir.c} triggering \option{-Wformat}
    is fixed.
  }
}
