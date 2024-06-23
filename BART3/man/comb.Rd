\name{comb}
\alias{comb}
\title{Generate a set of combinations \code{M} choose \code{K}.}

\description{ Generically, generate a set of unordered combinations
  without limits except for \eqn{K<=M}.
}

\usage{
comb(M, K, S=1:M)
}

\arguments{
\item{M}{The number of items. }
\item{K}{The number of items to choose.}
\item{S}{The items are represented by the integers \code{1:M} by
  default.  However, any arbitrary vector of integers length \code{M}
  will suffice.
}
}

\details{ A matrix of combinations is returned.  }

\seealso{
   \code{\link{SHAP}}
}
 
\keyword{marginal effects}

