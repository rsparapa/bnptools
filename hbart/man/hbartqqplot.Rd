\name{hbartqqplot}
\alias{hbartqqplot}
\docType{package}
\title{Predictive qqplot for heterbart}
\description{
Given results from predict.hbart, gets draws from the predictive distribution at each x
and then computes the empirical inverse cdf to get draws which would be uniform if the predicitive
were the true distribution.  Then draws a qqplot against the uniform.
In large enough samples, if the model is correct, the qqplot should look like like a straight line
with intercept 0 and slope 1.  In small samples, we expect the predictive to be more spread out than the true
distribution, even if the model is correct.
}
\usage{
hbartqqplot(y,rbmod,nunif=10000,linecolor="red",linewd=3,...)
}
\arguments{
\item{y}{y values corresponding to x values given to predict.hbart.}
\item{rbmod}{Output list from predict.hbart.}
\item{nunif}{Number of uniform(0,1) draws used in constructing the qqplot.}
\item{linecolor,linewd}{Line color and width for (0,1) line.}
\item{...}{Arguments passed on to stats::qqplot.}
}
\value{quantiles of y in draws from the predictive (conditional on each x value).}
\examples{
##################################################
## please see vignette and/or www.rob-mcculloch.org for more realistic examples
##################################################

## get simulated data
data(simdat)

##get hbart run on the simulated data
data(hbartonsimd)

##  Predictive quantile-quantile plot
temp = hbartqqplot(simdat$yp,hbartonsimd,xlab="predictive quantile",ylab="uniform",
                      cex.axis=1.4,cex.lab=1.2)
}

