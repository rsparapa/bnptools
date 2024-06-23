\name{read.trees}
\alias{read.trees}
\title{Reading in an ensemble of trees from BART}

\description{ BART encodes trees in ASCII strings for convenience.
  This function reads in the ensemble provided by BART posterior MCMC draws.
}

\usage{
read.trees(treedraws, ## treedraws item returned from BART
           x.train=matrix(nrow=0, ncol=0),
                      ## x.train to estimate coverage
           call=FALSE,## default to R vs. C++ code
           cutpoints=NULL,
           trees=NULL)
}

\arguments{
\item{treedraws}{The \code{treedraws} object from a BART fit object.
}
\item{x.train}{A training matrix.  Likely the same matrix that
  was used in the BART fit.  However, it may be from a different
  reference population that is similarly representative.}
\item{call}{Whether to call the R or C++ implementation.
  \code{FALSE} defaults to calling R.
}
\item{cutpoints}{A substitute set of cutpoints if provided.}
\item{trees}{A substitute ASCII string for the trees if provided.}
}

\details{ An object for the trees is returned that is appropriate for
the \code{EXPVALUE} function. }

\seealso{
   \code{\link{EXPVALUE}}
}
 
\keyword{marginal effects}

