\name{hbartonsimd}
\alias{hbartonsimd}
\docType{data}
\title{
hbart run on simulated data
}
\description{
predict.hbart results for simulated data.
}
\usage{data("hbartonsimd")}
\format{
  \describe{
    \item{\code{hbartonsimd}}{returned list from a call to predict.hbart on the simulated data.}
  }
}
\details{
The data hbartonsimd is the results of an hbart run on the simulated data in simdat.

The code for the hbart run is:

## load data \cr
data(simdat) \cr
attach(simdat) #some people think this is a bad idea \cr

## run hbart \cr
set.seed(99) \cr
res = hbart(x,y,nskip=500,ndpost=100,nadapt=500,adaptevery=50) \cr
hbartonsimd = predict(res,xp) #get prediction for test x in xp \cr
}
\examples{
## load simulated data and hbart run.
data(hbartonsimd)
data(simdat)

## plot data and x vs f(x), x vs f(x) +/- 2s(x), x in test data, true and estimated
## data
plot(simdat$xp,simdat$yp)
## true
lines(simdat$xp,simdat$fxp,col="blue",lty=2,lwd=2)
lines(simdat$xp,simdat$fxp+2*simdat$sxp,col="blue",lty=2,lwd=2)
lines(simdat$xp,simdat$fxp-2*simdat$sxp,col="blue",lty=2,lwd=2)
##estimated
mhat = hbartonsimd$mmean; shat = hbartonsimd$smean
lines(simdat$xp,mhat,col="red",lty=1,lwd=2)
lines(simdat$xp,mhat + 2*shat,col="red",lty=1,lwd=2)
lines(simdat$xp,mhat - 2*shat,col="red",lty=1,lwd=2)

## note that you can get "nicer" looking fits by 
## (i) running hbart longer (e.g. ndpost=500), 
## (ii) using numcut=1000,k=5 in hbart.
}
\keyword{datasets}
