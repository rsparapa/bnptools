\name{bartModelMatrix}
\alias{bartModelMatrix}
\title{Create a matrix out of a vector or data.frame}
\description{
  The compiled functions of this package operate on matrices in memory.  Therefore,
  if the user submits a vector or data.frame, then this function converts
  it to a matrix.  Also, it determines the number of cutpoints necessary
  for each column when asked to do so.
}
\usage{
bartModelMatrix(X, numcut=0L, usequants=FALSE, type=7, rm.const=FALSE,
                cont=FALSE, xicuts=NULL, rm.vars=NULL)
}

\arguments{
  \item{X}{A vector or data.frame to create the matrix from. }
  \item{numcut}{The maximum number of cutpoints to consider.
    If \code{numcut=0}, then just return a matrix; otherwise,
    return a list. }
  \item{usequants}{ If \code{usequants} is \code{FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, quantiles are used for the cutpoints. }
  \item{type}{ Determines which quantile algorithm is employed.}
  \item{rm.const}{ Whether or not to remove constant variables.}
  \item{cont}{ Whether or not to assume all variables are continuous.}
  \item{xicuts}{ To specify your own cut-points, use the \code{xicuts}
    argument.}
  \item{rm.vars}{ The variables that you want removed. }
}

\seealso{
   \code{\link{xicuts}}
 }

 \examples{

set.seed(99)

a <- rbinom(10, 4, 0.4)

table(a)

x <- runif(10)

df <- data.frame(a=factor(a), x=x)

b <- bartModelMatrix(df)

b

b <- bartModelMatrix(df, numcut=9)

b

b <- bartModelMatrix(df, numcut=9, usequants=TRUE)

b

\dontrun{
    f <- bartModelMatrix(as.character(a))
}
}
