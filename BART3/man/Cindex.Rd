\name{Cindex}
\alias{Cindex}
\title{C index for assessing survival analysis concordance}

\description{ The concordance probability, as estimated via the C index,
  is the most popular frequentist assessment of survival analyis model
  performance.  }

\usage{
Cindex(risk, times, delta=NULL)
}

\arguments{
\item{risk}{A vector of risk, i.e., higher (lower)
    risk is associated with shorter (longer) survival times.
}
\item{times}{A vector of the survival times.}
\item{delta}{An event indicator: 1, event; 0, censoring.  Defaults to
  all events.}
}

\details{ The C index is returned.  }

\references{
  FE Harrell Jr., KL Lee, RM Califf, DB Pryor and RA Rosatie (1984).
  Regression modelling strategies for improved prognostic prediction.
  \emph{Stat in Med}, \bold{3}, 143-152.
}

\seealso{
   \code{\link{surv.bart}}
}
 
\keyword{survival analysis}

