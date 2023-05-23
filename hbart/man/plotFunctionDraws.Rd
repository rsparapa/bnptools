\name{plotFunctionDraws}
\alias{plotFunctionDraws}
\docType{package}
\title{Plot matrix of function draws evaluated on a set of x}
\description{
Given draws of a function \eqn{f_d}{fd}, \eqn{d=1,2,\ldots,D}{d=1,2,...D} and a set of \eqn{x} vectors
\eqn{x_j, j=1,2,\ldots,J}{x(j), j=1,2,...J},  we have a \eqn{D \times J}{D x J} matrix of evaluations 
whose \eqn{(d,j)} element
is \eqn{f_d(x_j)}{fd(x(j))}, the \eqn{d^{th}}{d th} draw of the function evaluated at the \eqn{j^{th}}{j th} \eqn{x}. 
This function plots the draws by plotting estimates of \eqn{f(x_j)} versus intervals for \eqn{f(x_j)}. 
The estimates are the mean of the \eqn{j^{th}} column
and the intervals are two quantiles of the \eqn{j^{th}} column (e.g 5\% and 95\%).
}
\usage{
plotFunctionDraws(fd,complevel=mean(fd),probs=c(.025,.975),
      xlab="posterior mean of function",ylab="posterior intervals",
      intervalcol="green",linecol="red",
      pts=NA,ptscol="blue", ptspch=1, ptscex=1, ...)
}
\arguments{
\item{fd}{
\eqn{D \times J}{D times J} matrix whose \eqn{(d,j)}  element is the \eqn{d^th} function
draw evaluated at the \eqn{j^{th}} \eqn{x}.
}
\item{complevel}{
A horizontal line is drawn a complevel to compare the intervals to.
}
\item{probs}{
The two quantiles used to construct the intervals.
}
\item{xlab}{Label for x axis.}
\item{ylab}{Label for y axis.}
\item{intervalcol}{Color to draw the intervals with.}
\item{linecol}{Color to draw the comparizon horizontal line with.}
\item{pts}{Add \eqn{(x_j,pts_j)} to the plot.  For example pts could be fitted values from
and alternative model such as the linear model.}
\item{ptscol}{Color to draw the points pts with.}
\item{ptspch}{plot charactor to plot the points pts with.}
\item{ptscex}{cex to plot the points pts with.}
\item{...}{Arguments passed on to call to graphics::plot.}
}
\value{NULL}
\examples{
##################################################
## please see vignette and/or www.rob-mcculloch.org for more realistic examples
##################################################

## get simulated data
data(simdat)

##get hbart run on the simulated data
data(rbartonsimd)

## plot function (f and s) draws
shat = sqrt(mean((simdat$yp-rbartonsimd$mmean)^2)) #overall estimate of sigma
lmfit = lm(y~x,data.frame(x=simdat$x,y=simdat$y)) 
yhatlm = predict(lmfit,data.frame(x=simdat$xp)) #fits from a linear model

#Now we use plotFunctionDraws to look at mdraws (left panel) and sdraws (right panel).

## in the mean inference, you can see that the linear model seem unlikely
## in the variance inference, you can see that the posteriors of s(x) are far from a constant value
par(mfrow=c(1,2))

## look at mean inference
plotFunctionDraws(rbartonsimd$mdraws,complevel=mean(simdat$y), probs=c(.05,.95),
   xlab=expression(hat(f)(x)), pts=yhatlm, ptscol="black",
   cex.lab=1.2, cex.axis=1.4, main="intervals for f(x)")

##look at the standard deviation inference
plotFunctionDraws(rbartonsimd$sdraws, complevel=shat, xlab=expression(hat(s)(x)),
   intervalcol="magenta", linecol="blue",
   cex.lab=1.2, cex.axis=1.4, main="intervals for s(x)")

}


