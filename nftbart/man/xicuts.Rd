\name{xicuts}
\title{Specifying cut-points for the covariates}
\alias{xicuts}
\description{This function allows you to create a list that specifies
  the cut-points for the covariates.}
\usage{
xicuts(x.train, transposed=FALSE, numcut=100)
}
\arguments{
  \item{x.train}{The training matrix to derive cut-points from.}
  \item{transposed}{Whether or not the matrix has been tranposed yet.}
  \item{numcut}{The number of cut-points to create.}
}
\details{The cut-points are generated uniformly from min. to max., i.e.,
  the distribution of the data is ignored.
}
\value{
An object is returned of type \code{BARTcutinfo} which is essentially a list.
}
