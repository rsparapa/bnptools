\name{srstepwise}
\alias{srstepwise}
\title{Stepwise Variable Selection Procedure for survreg and lm}
\usage{
srstepwise(x, y, delta=1, sle = 0.15, sls = 0.15, dist='lognormal')
}
\arguments{
\item{x}{Matrix of variables to consider.}
  
\item{y}{If \code{dist=='normal'}, then this is a continuous outcome
  for a linear regression.  Otherwise, this is the time to an event
  for a parametric survival regression.}

\item{delta}{The event indicator: 1 for event, 0 for no event.}

\item{sle}{The chosen significance level for entering.}

\item{sls}{The chosen significance level for staying.}

\item{dist}{If \code{dist=='normal'}, then perform a linear regression
  with \code{lm}.
  Otherwise, this is the distribution to be used by \code{survreg}.}
}
\value{
Returns a list of indices of variables which have entered and stayed.
}
\description{
This stepwise variable selection procedure can be applied to obtain the
best candidates for a \code{lm} or \code{survreg} fit.
}
\details{
  Unfortunately, no stepwise procedure exists for \code{survreg} models.
  Therefore, we provide this brute force method and we provide the
  option to use the same function for \code{lm}.
}

\examples{

names. <- names(lung)[-(2:3)]
status1 <- ifelse(lung$status==2,1,0)
X <- as.matrix(lung)[ , names.]
vars=srstepwise(X, lung$time, status1)
print(names.[vars])

N = 1000
P = 5       
set.seed(12)
X=matrix(runif(N*P, -2, 2), N, P)
y = X[ , 1]^3
vars=srstepwise(X, y, dist='normal')
print(vars)

}
\seealso{
\link{lung}
}
