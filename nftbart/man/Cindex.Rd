\name{Cindex}
\alias{Cindex}
\title{Calculate the C-index for survival analysis.}

\description{ The C-index for survival analysis is the corollary of the
c statistic (the area under the Receiver Operating Characteristic curve)
for binary outcomes. As a probability of concordance, the higher is the
C-index, the better is the model discrimination vs. lesser probability
values. }

\usage{
Cindex(risk, times, delta=NULL)
}

\arguments{
  \item{risk}{A vector or prognostic risk scores. }
  \item{times}{A vector of failure times. }
  \item{delta}{ The corresponding failure time status code: 0,
    right-censored; 1, failure;
    or 2, left-censored. Defaults to all failures if not specified.}
}
\value{
The return value is the calculated C-index.
}
\references{
Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA. (1982)
  Evaluating the yield of medical tests. JAMA, May 14;247(18):2543-6.
}
\seealso{
   \code{\link{predict.nft}}
 }

\examples{

\dontrun{

}

}
