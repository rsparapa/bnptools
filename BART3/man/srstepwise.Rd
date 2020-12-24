\name{srstepwise}
\alias{srstepwise}
\title{Stepwise Variable Selection Procedure for survreg}
\usage{
srstepwise(x, y, delta=1, sle = 0.15, sls = 0.15, dist='lognormal')
}
\arguments{
\item{x}{Matrix of variables to consider.}
  
\item{y}{The outcome which is either continuous or time to an event.}

\item{delta}{The event indicator: 1 for an event, 0 for censoring.
Defaults to \code{1} for a continuous outcome with no censoring. }

\item{sle}{The chosen significance level for entering.}

\item{sls}{The chosen significance level for staying.}

\item{dist}{The distribution to be used by \code{survreg}.}
}
\value{
Returns a list of indices of variables which have entered and stayed.
}
\description{
This stepwise variable selection procedure can be applied to obtain the
best candidates for a \code{survreg} fit.
}
\details{
  Unfortunately, no stepwise procedure exists for \code{survreg} models.
  Therefore, we provide this brute force method.
}

\examples{

names. <- names(lung)[-(2:3)]
status1 <- ifelse(lung$status==2,1,0)
X <- as.matrix(lung)[ , names.]
vars=srstepwise(X, lung$time, status1)
print(names.[vars])

}
\seealso{
\link{lung}
}
