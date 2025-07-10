%% \section{Changes in version x.y.z (yyyy-mm-dd)}{
%%   \subsection{Significant User-Visible Changes}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Deprecated and Defunct}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Bug Fixes}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{New Features}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Installation}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Documentation}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Internals}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%% }

\name{NEWS}
\title{News for \R{} Package \pkg{fastbeta}}
\newcommand{\GH}{\href{https://github.com/davidearn/fastbeta/issues/#1}{GH##1}}

\section{Changes in version 0.4.1 (2025-07-10)}{
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item \code{sir.aoi(\dots, root = \var{function})} is a multiple
      time series object storing the solution until the last root or
      time point.  Attribute \code{root.info} stores details about the
      system at the roots in the format \code{list(tau, state)}.
      Previously, \code{sir.aoi} discarded the solution and returned
      only details about the system at the first root.

      \item \code{sir.aoi(\dots, root = \var{function})} no longer has
      an attribute \code{curvature} storing the curvature of \code{Y} at
      the root.
    }
  }
  \subsection{Deprecated and Defunct}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item The expression for \eqn{\text{d}Y/\text{d}\tau}{dY/dtau}
      used by \code{sir.aoi} was wrong.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item \code{sir.aoi(\dots, aggregate = TRUE)} has additional
      columns storing the force of infection, incidence, and the
      curvature of \code{Y}.

      \item \code{sir.aoi} accepts new argument \code{eps} in order to
      support equations with equal, nonzero rates of birth and death.

      \item \code{sir.aoi} accepts new arguments \code{F} and
      \code{Fargs}, enabling the user to set a forcing function.

      \item \code{sir.aoi} accepts new arguments \code{H} and
      \code{Hargs}, enabling the user to set a susceptible
      hetereogeneity function.

      \item \code{sir.aoi} accepts new arguments \code{root.max} and
      \code{root.break}, enabling the user to set the number of roots
      sought and indicate whether the solver should stop once that
      number of roots is reached.  The default behaviour is to stop once
      one root is found and agrees with the previous behaviour.

      \item The method for generic function \code{summary} and class
      \code{sir.aoi} accepts new argument \code{name}, enabling
      approximation of tail exponents of different variables.

      \item Function \code{cbind.ts} for working around \PR{18583}.

      \item New data set \code{influenza}, used in an example of
      \code{deconvolve} usage.
    }
  }
  \subsection{Installation}{
    \itemize{
      \item .
    }
  }
  \subsection{Documentation}{
    \itemize{
      \item .
    }
  }
  \subsection{Internals}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.4.0 (2025-04-28)}{
  \subsection{Deprecated and Defunct}{
    \itemize{
      \item \code{seir.canonical} is removed in favour of more general
      \code{sir.aoi}.  The method for generic function \code{summary}
      and class \code{seir.canonical} is removed in favour of a method
      for class \code{sir.aoi}.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Function \code{sir.aoi} generalizing \code{seir.canonical}
      in order to support rates of transmission and recovery structured
      by age of infection.
    }
  }
}

\section{Changes in version 0.3.2 (2025-03-23)}{
  \subsection{New Features}{
    \itemize{
      \item Function \code{seir.canonical} providing a nondimensional
      interface to a special case (constant transmission, without vital
      dynamics, without loss of immunity) of \code{seir} usage.  The
      nondimensional parameters are the basic reproduction number and
      the ratio of the mean latent period and mean generation interval.

      \item A method for generic function \code{summary} and class
      \code{seir.canonical} approximating the tail exponents of
      prevalence.
    }
  }
}

\section{Changes in version 0.3.1 (2024-11-22)}{
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item Formal argument \code{init} of \code{fastbeta},
      \code{fastbeta.bootstrap}, \code{ptpi}, and \code{seir} is
      repositioned.

      \item Formal arguments \code{sigma} and \code{m} of
      \code{fastbeta}, \code{fastbeta.bootstrap}, \code{ptpi}, and
      \code{seir} have constant default values (\code{sigma = 1},
      \code{m = 1L}) not depending on the values of other arguments.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Entry \code{[2, 2]} of the Jacobian matrix constructed by
      \code{seir(stochastic = FALSE, useCompiled = FALSE)} was wrong.

      \item \code{seir} handles more comprehensively the situation in
      which the number of rows of the matrix returned by \pkg{deSolve}
      function \code{lsoda} is less than \code{length.out}.  It now
      distinguishes termination due to found roots from termination due
      to an unsuccessful solver call.

      \item Warnings caught by \option{-Wconversion} are squashed.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item \code{seir} gains logical argument \code{aggregate}
      indicating whether latent and infectious compartments should be
      aggregated.

      \item \code{seir} arguments \code{nu} and \code{mu} gain
      default values.  In both cases, it is \code{function (t) 0},
      indicating no vital dynamics.

      \item \code{seir.R0}, \code{seir.ee}, and \code{seir.jacobian} for
      calculating the basic reproduction number, endemic equilibrium,
      and Jacobian matrix.

      \item \code{fastbeta.matrix} for calculating lower triangular
      coefficient matrix.
    }
  }
  \subsection{Installation}{
    \itemize{
      \item Compilation now uses \option{-DSTRICT_R_HEADERS} and
      \option{-DR_NO_REMAP}.
    }
  }
}

\section{Changes in version 0.3.0 (2024-05-17)}{
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item Argument \code{constants} of functions \code{sir},
      \code{fastbeta}, \code{fastbeta.bootstrap}, and \code{ptpi} is
      replaced by arguments \code{sigma}, \code{gamma}, \code{delta},
      \code{init}, \code{m}, and \code{n}.

      \item \code{sir} is now a wrapper calling new \code{seir}.

      \item Arguments \code{a} and \code{b} of function \code{ptpi}
      are replaced by arguments \code{start} and \code{end} so that
      the interface is more similar to that of \code{window}.

      \item \code{ptpi} returns \code{list(value, diff, iter, x)},
      no longer \code{list(value, delta, iter, X)}.

      \item Simulations \code{sir.e01} and \code{sir.e02} are replaced
      by \code{seir.ts01} and \code{seir.ts02}.  The replacements have
      a latent compartment and a different set of attributes reflecting
      the API of \code{seir}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item More simulations set \code{tl.params = list(epsilon = .)}
      to mitigate noise attributable to too big leaps in the adaptive
      tau-leaping algorithm.

      \item \code{plot} method for class \code{fastbeta.bootstrap}
      uses \code{dev.hold} and \code{dev.flush} from \pkg{grDevices}
      to avoid \dQuote{incremental} graphical output on some devices.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Compartmental model is generalized to support multiple
      infected compartments, \eqn{m} latent and \eqn{n} infectious,
      resulting in several backwards incompatible changes to the API;
      see below.
    }
  }
  \subsection{Installation}{
    \itemize{
      \item \pkg{adaptivetau} and \pkg{deSolve} are now suggested,
      not imported, as they are needed only for simulation.
      Usage of simulated data sets is now everywhere conditional on
      successful loading of the dependency.
    }
  }
}

\section{Changes in version 0.2.0 (2023-12-11)}{
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item \code{sir}
      expects \code{constants = c(S0, I0, R0, gamma, delta)},
      no longer \code{constants = c(gamma, S0, I0, R0)}.

      \item \code{fastbeta} and \code{fastbeta.bootstrap}
      expect \code{constants = c(S0, I0, R0, gamma, delta)},
      no longer \code{constants = c(gamma, S0, I0)}.

      \item \code{ptpi} gains argument \code{constants} and loses
      argument \code{start}.
      It expects \code{constants = c(Sa, Ia, Ra, gamma, delta)},
      no longer \code{start = Sa}.

      \item \code{fastbeta} returns (modulo attributes)
      \code{cbind(S, I, R, beta)}, no longer \code{cbind(beta, S, I)}.

      \item \code{ptpi(\dots)[[1L]]} gives the state at time \code{0}
      or at time \code{a} conditional on new argument \code{backcalc}.
      The default is \code{a}, no longer (implicitly) \code{0}.

      \item \code{ptpi(\dots)[[4L]]} has dimensions
      \code{c(b-a+1, 3, iter)}, no longer \code{c(b-a+1, iter)}.

      \item Simulation \code{sir.e01} is generated by the same
      parameters but is no longer stochastic.
      Simulation \code{sir.e02} gives the stochastic variant.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Compartmental model is generalized to support loss of
      immunity at rate \eqn{delta}, resulting in several backwards
      incompatible changes to the API; see below.
    }
  }
  \subsection{Installation}{
    \itemize{
      \item \file{data/*.R} are no longer copied to \file{inst/scripts}
      at install time by a \command{make} rule in \file{src/Makevars}.
      There is no need: the scripts are preserved due to the setting of
      \code{BuildResaveData: no}.  Hence \file{src/Makevars} is removed.
    }
  }
}

\section{Changes in version 0.1.1 (2023-11-28)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item Typo in \file{src/sir.c} triggering \option{-Wformat}
      is fixed.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Changes are logged in \file{inst/NEWS.Rd}.

      \item \file{DESCRIPTION} sets \code{BuildResaveData: no},
      so that \command{R CMD build} works without
      \option{--no-resave-data} when \pkg{fastbeta} is not installed.
    }
  }
  \subsection{Installation}{
    \itemize{
      \item \file{data/*.R} are copied into \file{inst/scripts}
      at install time by a \command{make} rule in \file{src/Makevars}.

      \item \code{Depends: R (>= 4.3.0)} for \verb{amsmath} in PDF
      output, which was not supported until 4.2.2.
    }
  }
}
