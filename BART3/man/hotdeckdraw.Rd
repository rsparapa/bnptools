\name{hotdeckdraw}
\title{Hot-deck missing imputation}
\alias{hotdeckdraw}
\description{This function imputes missing data with hot-decking.}
\usage{
hotdeckdraw(x.train, y.train, x.test=matrix(0, 0, 0), impute.mult=NULL)
}
\arguments{
  \item{x.train}{The training matrix.}
  \item{y.train}{The outcome to measure nearness.}
  \item{x.test}{The testing matrix, if given.}
  \item{impute.mult}{An index of the columns to avoid
    imputing which will be handled by BART internally.}
}
\details{
  Hot-decking was a method commonly employed with US Census data
  in the early computing era.  For a particular respondent, missing
  data was imputed by randomly selecting from the responses of their
  neighbors since it is assumed that the values are likely similar.
  In our case, we measure nearness via the outcome \code{y.train}.
  However, this can only be done for training since testing may,
  or may not, be available at the time of fitting.  If the training and testing
  matrices are the same, then they receive hot-deck imputation;
  otherwise testing uses cold-decking.
}
\value{
  \item{x.train}{The imputed training matrix.}
  \item{x.test}{The imputed testing matrix.}
  \item{miss.train}{A summary of the missing variables for training.}
  \item{miss.test}{A summary of the missing variables for testing.}
  \item{impute.flag}{Whether \code{impute.mult} columns were, or were
    not, imputed.}
  \item{same}{Whether \code{x.train} and \code{x.test} are identical.}
}
