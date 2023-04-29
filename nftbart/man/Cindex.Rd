\name{Cindex}
\alias{Cindex}
\alias{concordance}
\title{Calculate the C-index/concordance for survival analysis.}

\description{ The C-index for survival analysis is the corollary of the
c statistic (the area under the Receiver Operating Characteristic curve)
for binary outcomes. As a probability, the higher is the C-index, the
better is the model discrimination vs. lesser probability
values. Similarly, the concordance is calculated like the C-index from
z-draws via the posterior predictive distribution restricted to the
horizon of the data (a la restricted mean survival time).}

\usage{
Cindex(risk, times, delta=NULL)

concordance(draws, times, delta=NULL)
}

\arguments{
  \item{risk}{A vector or prognostic risk scores. }
  \item{draws}{A vector of draws via
the posterior predictive distribution restricted to the horizon of
the data (a la restricted mean survival time). }
  \item{times}{A vector of failure times. }
  \item{delta}{ The corresponding failure time status code: 0,
    right-censored; 1, failure;
    or 2, left-censored. Defaults to all failures if not specified.}
}
\value{
The return value is the calculated C-index/concordance.
}
\references{
Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA. (1982)
  Evaluating the yield of medical tests. JAMA, May 14;247(18):2543-6.
}
\seealso{
   \code{\link{predict.nft}}
 }

\examples{

data(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1

## this study reports time in days
times=lung$time
times=times/7  ## weeks

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])
## lung$sex:        Male=1 Female=2

\dontrun{
    set.seed(99)
    post=nft(x.train, times, delta, K=0)
    pred=predict(post, x.train, XPtr=TRUE, seed=21)
    print(Cindex(pred$logt.test.mean, times, delta))
}

}
