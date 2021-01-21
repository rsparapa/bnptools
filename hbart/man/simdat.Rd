\name{simdat}
\alias{simdat}
\docType{data}
\title{
Simulated Example
}
\description{
Simulated data with nonlinear mean and heteroskedasticity.
}
\usage{data("simdat")}
\format{
  \describe{
    \item{\code{x}}{simulated train x values}
    \item{\code{y}}{simulated train y values}
    \item{\code{xp}}{simulated test xp values}
    \item{\code{yp}}{simulated test yp values}
    \item{\code{fx}}{true f evaluated on train x}
    \item{\code{sx}}{true s evaluated on train x}
    \item{\code{fxp}}{true f evaluated on test xp}
    \item{\code{sxp}}{true s evaluated on test xp}
  }
}
\details{
The simulated data in simdat was generated using the code: 

##simulate data  \cr
set.seed(99)

# train data \cr
n=500 #train data sample size  \cr
p=1 #just one x \cr
x = matrix(sort(runif(n*p)),ncol=p) #iid uniform x values \cr
fx = 4*(x[,1]^2) #quadratric function f \cr
sx = .2*exp(2*x[,1]) # exponential function s \cr
y = fx + sx*rnorm(n) # y = f(x) + s(x) Z \cr

#test data (the p added to the variable names is for predict) \cr
np=1000 #test data sample size  \cr
xp = matrix(sort(runif(np*p)),ncol=p) \cr
fxp = 4*(xp[,1]^2) \cr
sxp = .2*exp(2*xp[,1]) \cr
yp = fxp + sxp*rnorm(np) \cr
}
\examples{
data(simdat)

## plot x vs y with f(x) and f(x) +/- 2s(x) for train and test simulated data
##train
plot(simdat$x,simdat$y,xlab="x",ylab="y")
##test
points(simdat$xp,simdat$yp,col="red",pch=2)
lines(simdat$xp,simdat$fxp,col="blue",lwd=2)
lines(simdat$xp,simdat$fxp+2*simdat$sxp,col="blue",lwd=2,lty=2)
lines(simdat$xp,simdat$fxp-2*simdat$sxp,col="blue",lwd=2,lty=2)
legend("topleft",legend=c("train","test"),pch=c(1,2),col=c("black","red"))
}
\keyword{datasets}
