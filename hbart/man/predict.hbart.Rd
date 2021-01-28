\name{predict.hbart}
\alias{predict.hbart}
%\docType{package}
\title{Drawing Posterior Predictive Realizations for hbart models.}
\description{
  The function \code{predict.hbart()} is the main function for drawing posterior predictive realizations at new inputs using a fitted model stored in a \code{hbart} object returned from \code{hbart()}.
}
\usage{
\method{predict}{hbart}(
object,
x.test=object$x.train,
tc=1,
fmean=mean(object$y.train),
q.lower=0.025,
q.upper=0.975,...)
}
\arguments{
  \item{object}{Object of type \code{hbart} from a previous call to \code{hbart()}}
  \item{x.test}{New input settings in the form of an \code{npred x p} matrix at which to construct predictions.  Defaults to the training inputs.}
  \item{tc}{Number of OpenMP threads to use for parallel computing.}
  \item{fmean}{Mean-centering vector for the training data.  Defaults to the value used by \code{hbart()} when fitting the model.  Usually should be left to the default.}
  \item{q.lower}{Lower quantile to return.}
  \item{q.upper}{Upper quantile to return.}
   \item{...}{not used.}
}
\details{
  \code{predict.hbart()} is the main function for calculating posterior predictions and uncertainties once a model has been fit by \code{hbart()}.  

  Returns an object of type \code{hbart} with the following entries.
}
\value{
  \item{mdraws}{Posterior realizations of the mean function, \eqn{f({\bf x})}{f(x)} stored in an \code{ndpost x npred} matrix,
    where ndpost is the number of kept MCMC draws in the \code{hbart} run.}
  \item{sdraws}{Posterior realizations of the standard deviation function, \eqn{s({\bf x})}{s(x)} stored in an \code{ndpost x npred} matrix,
    where ndpost is the number of kept MCMC draws in the \code{hbart} run.}
  \item{mmean}{Posterior predictive mean of \eqn{f({\bf x})}{f(x)}.}
  \item{smean}{Posterior predictive mean of the standard deviation, \eqn{s({\bf x})}{s(x)}.}
  \item{msd}{Posterior standard deviation of the mean, \eqn{f({\bf x})}{f(x)}.}
  \item{ssd}{Posterior standard deviation of the standard devation, \eqn{s({\bf x})}{s(x)}.}
  \item{m.5}{Posterior median of the mean function realizations, \eqn{f({\bf x})}{f(x)}.}
  \item{m.lower}{Posterior \code{q.lower} quantile of the mean function realizations.}
  \item{m.upper}{Posterior \code{q.upper} quantile of the mean function realizations.}
  \item{s.5}{Posterior median of the standard deviation function realizations, \eqn{s({\bf x})}{s(x)}.}
  \item{s.lower}{Posterior \code{q.lower} quantile of the standard deviation function realizations.}
  \item{s.upper}{Posterior \code{q.upper} quantile of the standard deviation function realizations.}
  \item{q.lower}{Lower quantile used in constructing the above.}
  \item{q.upper}{Upper quantile used in constructing the above.}
}
\references{
  Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (1998)
  Bayesian CART model search. 
  \emph{Journal of the American Statistical Association}, \bold{93}, 935--948.

  Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (2010)
  BART: Bayesian additive regression trees. 
  \emph{The Annals of Applied Statistics}, \bold{4}, 266--298.

  Pratola, Matthew T. (2016)
  Efficient Metropolis Hastings proposal mechanisms for Bayesian regression tree models.
  \emph{Bayesian analysis}, \bold{11}, 885--911.
 
  Pratola, Matthew T., Chipman, Hugh A., George, Edward I. and McCulloch, Robert E. (2017)
  Heteroscedastic BART Using Multiplicative Regression Trees.
  \emph{arXiv preprint}, \bold{arXiv:1709.07542}, 1--20.
}
\author{
Matthew T. Pratola <mpratola@stat.osu.edu> [aut, cre, cph], Robert E. McCulloch <robert.e.mcculloch@gmail.com> [aut, cre, cph], Hugh Chipman <hugh.chipman@gmail.com> [aut, cph]
Maintainer: Matthew T. Pratola <mpratola@stat.osu.edu>, Robert E. McCulloch <robert.e.mcculloch@gmail.com>
}
\seealso{
\code{\link{hbart}}
}
\examples{
##################################################
## This is just a stub (runs fast) example for testing.
##  For more realistic examples, please see:
##   (i) the vignette at www.rob-mcculloch.org
##   (ii) the example simulated data (see ?simdat)
##        and the longer run in ?hbartonsimd, 
##        where a saved run of hbart is run on simdat is plotted. 
##################################################

##simulate data
set.seed(99)

# train data
n=500 #train data sample size
p=1 #just one x
x = matrix(sort(runif(n*p)),ncol=p) #iid uniform x values
fx = 4*(x[,1]^2) #quadratric function f
sx = .2*exp(2*x[,1]) # exponential function s
y = fx + sx*rnorm(n) # y = f(x) + s(x) Z

#test data (the p added to the variable names is for predict)
np=500 #test data sample size
xp = matrix(sort(runif(np*p)),ncol=p)
fxp = 4*(xp[,1]^2)
sxp = .2*exp(2*xp[,1])
yp = fxp + sxp*rnorm(np)

##run hbart MCMC
# The number of interations is kept small to make example run,
##!!!!  REAL APPLICATIONS MAY NEED LONGER RUNS !!!!
#   nskip: burn in draws,
#   ndpost:kept draws,
#   nadapt: initial draws to tune MCMC,
#   numcut: number of cutpoints used for each x
#   k: bigger k gives smoother f (default is 2)
set.seed(19)
res = hbart(x,y,nskip=10,ndpost=20,nadapt=0,numcut=1000,k=5) #again, this is way too short a run!!!
## now predict to get inference
resp = predict(res,x.test=xp)

##check out of sample fit
cat("out of sample cor(f,fhat) is ",cor(fxp,resp$mmean),"\n")
cat("out of sample cor(s,shat) is ",cor(sxp,resp$smean),"\n")

##plot estimated vs. true
##plot the data
plot(xp,yp,cex.axis=1.5,cex.lab=1.5)
lines(xp,fxp,col="blue")
lines(xp,fx+2*sxp,col="blue",lty=2)
lines(xp,fxp-2*sxp,col="blue",lty=2)

## add the fit
lines(xp,resp$mmean) #estimate of f
lines(xp,resp$mmean+2*resp$smean) #estimate of sd
lines(xp,resp$mmean-2*resp$smean) #estimate of sd
}
