\name{nft}
\alias{nft}
%\docType{package}
\title{Fitting NFT BART models.}
\description{
  The function \code{nft()} is the main function for fitting Bayesian Regression Tree models, including single-tree models, Bayesian Additive Regression Tree (BART) models and Heteroscedastic BART models. 
}
\usage{
nft(x.train, y.train, x.test=matrix(0.0,0,0), ntree=200, ntreeh=40, ndpost=1000, 
     nskip=100, k=2, power=2.0, base=.95, tc=1, sigmav=rep(1,length(y.train)), 
     fmean=mean(y.train), overallsd = sd(y.train), overallnu=10, 
     chv = cor(x.train,method="spearman"), pbd=.7, pb=.5, stepwpert=.1, probchv=.1,
     minnumbot=5, printevery=100, numcut=100, xicuts=NULL, nadapt=1000, adaptevery=100,
     summarystats=FALSE)
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
  \code{nft()} is the main model fitting function for continuous response data.  The most general form of the model allowed is
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
  \item{res}{Fitted model object of S3 class nft.}
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
\code{\link{nft}},\code{\link{predict.nft}}
}
\examples{

\dontrun{

B=8
##B=getOption('mc.cores', 1)
##figures = getOption('figures', default='NONE')

## load survival package for the advanced lung cancer example
data(lung)

N=length(lung$status)

## if physician's KPS unavailable, then use the patient's
h=which(is.na(lung$ph.karno))
lung$ph.karno[h]=lung$pat.karno[h]

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1
table(delta)

## this study reports time in days rather than weeks or months
times=lung$time
times=times/7  ## weeks
summary(times)

## matrix of observed covariates
x.train=cbind(lung$sex, lung$age, lung$ph.karno)
## lung$sex:        Male=1 Female=2
## lung$age:        Age in years
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician
dimnames(x.train)[[2]]=c('M(1):F(2)', 'age(39:82)', 'ph.karno(50:100:10)')

file.='lung.rds'
if(file.exists(file.)) {
    post=readRDS(file.)
    XPtr=FALSE
} else {
    set.seed(99)
    post=nft(x.train, times, delta, tc=B, K=0)
    XPtr=TRUE
    saveRDS(post, file.)
}

x.test = rbind(x.train, x.train)
x.test[ , 1]=rep(1:2, each=N)
events=seq(5, 150, 5)
K=length(events)

pred = predict(post, x.test, K=K, events=events,
               XPtr=XPtr, tc=B, FPD=TRUE)

plot(c(0, events), c(1, pred$surv.fpd.mean[1:K]), type='l', col=4,
     ylim=0:1, 
     xlab=expression(italic(t)),
     ylab=expression(italic(S)(italic(t), italic(x))))
lines(c(0, events), c(1, pred$surv.fpd.upper[1:K]), lty=2, lwd=2, col=4)
lines(c(0, events), c(1, pred$surv.fpd.lower[1:K]), lty=2, lwd=2, col=4)
lines(c(0, events), c(1, pred$surv.fpd.mean[K+1:K]), lwd=2, col=2)
lines(c(0, events), c(1, pred$surv.fpd.upper[K+1:K]), lty=2, lwd=2, col=2)
lines(c(0, events), c(1, pred$surv.fpd.lower[K+1:K]), lty=2, lwd=2, col=2)
legend('topright', c('Adv. lung cancer\nmortality example',
                     'M', 'F'), lwd=2, col=c(0, 4, 2), lty=1)
##dev.copy2pdf(file='lung.pdf')


}
}
