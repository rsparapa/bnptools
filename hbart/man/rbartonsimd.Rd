\name{rbartonsimd}
\alias{rbartonsimd}
\docType{data}
\title{
rbart run on simulated data
}
\description{
predict.rbart results for simulated data.
}
\usage{data("rbartonsimd")}
\format{
  \describe{
    \item{\code{rbartonsimd}}{returned list from a call to predict.rbart on the simulated data.}
  }
}
\details{
The data rbartonsimd is the results of an rbart run on the simulated data in simdat.

The code for the rbart run is:

## load data \cr
data(simdat) \cr
attach(simdat) #some people think this is a bad idea \cr

## run rbart \cr
set.seed(99) \cr
res = rbart(x,y,nskip=500,ndpost=100,nadapt=500,adaptevery=50) \cr
rbartonsimd = predict(res,xp) #get prediction for test x in xp \cr
}
\examples{
## load simulated data and rbart run.
data(rbartonsimd)
data(simdat)

## plot data and x vs f(x), x vs f(x) +/- 2s(x), x in test data, true and estimated
## data
plot(simdat$xp,simdat$yp)
## true
lines(simdat$xp,simdat$fxp,col="blue",lty=2,lwd=2)
lines(simdat$xp,simdat$fxp+2*simdat$sxp,col="blue",lty=2,lwd=2)
lines(simdat$xp,simdat$fxp-2*simdat$sxp,col="blue",lty=2,lwd=2)
##estimated
mhat = rbartonsimd$mmean; shat = rbartonsimd$smean
lines(simdat$xp,mhat,col="red",lty=1,lwd=2)
lines(simdat$xp,mhat + 2*shat,col="red",lty=1,lwd=2)
lines(simdat$xp,mhat - 2*shat,col="red",lty=1,lwd=2)

## note that you can get "nicer" looking fits by 
## (i) running rbart longer (e.g. ndpost=500), 
## (ii) using numcut=1000,k=5 in rbart.
}
\keyword{datasets}
