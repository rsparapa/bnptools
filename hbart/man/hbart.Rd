\name{hbart}
\alias{hbart}
%\docType{package}
\title{Fitting Bayesian Regression Tree models supported by hbart.}
\description{
  The function \code{hbart()} is the main function for fitting Bayesian Regression Tree models, including single-tree models, Bayesian Additive Regression Tree (BART) models and Heteroscedastic BART models.  \code{hbart()} maintains some degree of backwards compatibility with \code{BayesTree::bart()}, while offering many new options.
}
\usage{
hbart(x.train, y.train, x.test=matrix(0.0,0,0), ntree=200, ntreeh=40,
      ndpost=1000, nskip=100, k=2, power=2.0, base=.95,
      tc=getOption('mc.cores', 1L), sigmav=rep(1,length(y.train)), 
      fmean=mean(y.train), overallsd = sd(y.train), overallnu=10, 
      chv = cor(x.train,method="spearman"), pbd=.7, pb=.5, stepwpert=.1,
      probchv=.1, minnumbot=5, printevery=100, numcut=100, xicuts=NULL,
      nadapt=1000, adaptevery=100, summarystats=FALSE)
}
\arguments{
  \item{x.train}{nxp matrix of predictor variables for the training data.}
  \item{y.train}{nx1 vector of the observed response for the training data.}
  \item{x.test}{mxp matrix of predictor variables for the test set. Deprecated.}
  \item{ntree}{Number of trees used for the mean model.}
  \item{ntreeh}{Number of trees used for the variance model.}
  \item{ndpost}{Number of iterations to run the MCMC algorithm after burn-in.}
  \item{nskip}{Number of MCMC iterations treated as burn-in and discarded.}
  \item{k}{Prior hyperparameter for the mean model.}
  \item{power}{Power parameter in the tree depth penalizing prior.}
  \item{base}{Base parameter in the tree depth penalizing prior.}
  \item{tc}{Number of OpenMP threads to use.}
  \item{sigmav}{Initialization of square-root of variance parameter.}
  \item{fmean}{Overall mean of the data for pre-centering the data before running the model.}
  \item{overallsd}{A rudimentary estimate of the process standard deviation. Used in calibrating the variance prior.}
  \item{overallnu}{Shape parameter for the variance prior.}
  \item{chv}{Predictor correlation matrix used as a pre-conditioner for MCMC change-of-variable proposals.}
  \item{pbd}{Probability of performing a birth/death proposal, otherwise perform a rotate proopsal.}
  \item{pb}{Probability of performing a birth proposal given that we choose to perform a birth/death proposal.}
  \item{stepwpert}{Initial width of proposal distribution for peturbing cutpoints.}
  \item{probchv}{Probability of performing a change-of-variable proposal.  Otherwise, only do a perturb proposal}.
  \item{minnumbot}{Minimum number of observations required in bottom (terminal) nodes.}
  \item{printevery}{Outputs MCMC algorithm status every printevery iterations.}
  \item{numcut}{Number of cutpoints to use for each predictor variable.}
  \item{xicuts}{More detailed construction of cutpoints can be specified using makecuts() and passed as an argument here.}
  \item{nadapt}{Number of MCMC iterations allowed for adaptive MCMC. These are also discarded.}
  \item{adaptevery}{Adapt MCMC proposal distributions every adaptevery iterations until the algorithm has run for nadapt iterations.}
  \item{summarystats}{Return detailed summary statistics about the fitting procedure.}
}
\details{
  \code{hbart()} is the main model fitting function for continuous response data.  The most general form of the model allowed is
  \eqn{Y({\bf x})=f({\bf x})+s({\bf x})Z}{Y(x)=f(x)+s(x)Z}
  where \eqn{Z} is \eqn{N(0,1)} and 
  \eqn{f({\bf x})=\sum_{j=1}^m g({\bf x};T_j,M_j)}{f(x)=sum g(x;T_j,M_j)} and
  \eqn{s({\bf x})=\prod_{j=1}^{m^\prime} h({\bf x};T^\prime_j,M^\prime_j)}{s(x)=prod g(x;T'_j,M'_j)},
  where the \eqn{g(\cdot;T_j,M_j)}{g(.;T_j,M_j)} represent additive tree components used for modeling the mean and
  \eqn{h(\cdot;T^\prime_j,M^\prime_j)}{h(.;T'_j,M'_j)} represent multiplicative tree components used for modeling the variance.

  The most common models to fit are a homoscedastic single-tree model, a homoscedastic BART model and a heteroscedastic BART model.

  For a BART model, set \code{pbd=c(0.7,0.0)} and \code{ntreeh=1}.  This forces a scalar (homoscedastic) variance term.

  For a single-tree model, set \code{pbd=(0.7,0.0)}, \code{ntreeh=1} and \code{ntree=1}. This forces the mean component to be modeled using only one tree.

  The heteroscedastic BART model is the default.

}
\value{
  \item{res}{Fitted model object of S3 class hbart.}
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
\code{\link{predict.hbart}}
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
