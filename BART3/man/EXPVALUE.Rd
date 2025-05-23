\name{EXPVALUE}
\alias{EXPVALUE}
\title{EXPVALUE function for BART}

\description{ The EXPVALUE function was originally proposed for a single
tree as a marginal expectation of \eqn{E(g(x.S, x.C; T, M)|x.S)} over
the training cohort where \eqn{S} are the covariate subset of interest
while \eqn{C} is the complement, i.e., \eqn{x.S} are settings of
interest and \eqn{x.C} are arbitrary.  Here we adapt it to an ensemble
of trees for BART: \eqn{E(f(x.S, x.C)|x.S)}. }

\usage{
EXPVALUE(trees,
         x.test,
         S,
         x.train)
}

\arguments{
\item{trees}{The forest of trees as decoded by \code{read.trees()}.
}

\item{x.test}{A matrix of settings.}
\item{S}{The indices of the columns of \code{x.test}
  that correspond to \code{x.train}, i.e., the
  first item of \code{S} is for the first column
  of \code{x.test} corresponding to the \code{S[1]}th column of
  \code{x.train}.}
\item{x.train}{A matrix of training values.}
}

\details{ The expectation is returned.  }

\references{
  SM Lundberg, G Erion, H Chen, A DeGrave, JM Prutkin, B Nair and R Katz
  (2020).
  From local explanations to global understanding with explainable AI
  for trees.
  \emph{Nature Machine Intelligence}, \bold{2}, 56-67.
}

\seealso{
   \code{\link{SHAP}}
}
 
\keyword{marginal effects}

