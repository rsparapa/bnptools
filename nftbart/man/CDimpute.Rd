\name{CDimpute}
\title{Cold-deck missing imputation}
\alias{CDimpute}
\description{This function imputes missing data.}
\usage{
CDimpute(x.train, x.test=matrix(0, 0, 0), impute.bin=NULL)
}
\arguments{
  \item{x.train}{The training matrix.}
  \item{x.test}{The testing matrix, if given.}
  \item{impute.bin}{An index of the columns to avoid
    imputing which will be handled by BART internally.}
}
\details{We call this method cold-decking in analogy to hot-decking.
  Hot-decking was a method commonly employed with US Census data
  in the early computing era.  For a particular respondent, missing
  data was imputed by randomly selecting from the responses of their
  neighbors since it is assumed that the values are likely similar.
  In our case, we make no assumptions about which values may, or may
  not, be nearby.  We simply take a random sample from the matrix
  rows to impute the missing data.  If the training and testing
  matrices are the same, then they receive the same imputation.
}
\value{
  \item{x.train}{The imputed training matrix.}
  \item{x.test}{The imputed testing matrix.}
  \item{miss.train}{A summary of the missing variables for training.}
  \item{miss.test}{A summary of the missing variables for testing.}
  \item{impute.flag}{Whether \code{impute.bin} columns were, or were
    not, imputed.}
  \item{same}{Whether \code{x.train} and \code{x.test} are identical.}
}
