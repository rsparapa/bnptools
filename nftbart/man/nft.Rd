\name{nft}
\alias{nft}
\title{Fit NFT BART models.}
\description{
  The \code{nft()} function is for fitting
  NFT BART (Nonparametric Failure Time
  Bayesian Additive Regression Tree) models. 
}
\usage{
nft(
    ## data
    x.train, times, delta=NULL, x.test=matrix(nrow=0, ncol=0),
    impute.bin=NULL, impute.prob=NULL,
    ## multi-threading
    tc=1, 
    ##MCMC
    nskip=1000, ndpost=2000, nadapt=1000, adaptevery=100,
    chv = cor(x.train, method="spearman"),
    pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
    stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
    minnumbot=c(5, 5),
    ## BART and HBART prior parameters
    ntree=c(50, 10), numcut=100, xicuts=NULL,
    power=c(2, 2), base=c(0.95, 0.95),
    ## f function
    k=5, sigmaf=NA, dist='weibull', 
    ## s function
    sigmav=NULL, overalllambda=NA, overallnu=10,
    ## survival analysis 
    K=100, events=NULL, 
    ## DPM LIO
    drawMuTau=1L, 
    alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
    neal.m=2, constrain=1, k0.draw=1, b0.draw=1,
    ## misc
    printevery=100
)

}
\arguments{
  \item{x.train}{nxp matrix of predictor variables for the training data.}
  \item{times}{nx1 vector of the observed times for the training data.}
  \item{delta}{nx1 vector of the time type for the training data:
  0, for right-censoring; 1, for an event; and, 2, for left-censoring.}
\item{x.test}{mxp matrix of predictor variables for the test set.}
  \item{impute.bin}{Indices of the columns of \code{x.train} to be imputed.}
  \item{impute.prob}{nx1 vector of prior probabilities for imputation.}
  \item{tc}{Number of OpenMP threads to use.}
  \item{nskip}{Number of MCMC iterations to burn-in and discard.}
  \item{ndpost}{Number of MCMC iterations kept after burn-in.}
  \item{nadapt}{Number of MCMC iterations for adaptation prior to burn-in.}
  \item{adaptevery}{Adapt MCMC proposal distributions every \code{adaptevery} iteration.}
    \item{chv}{Predictor correlation matrix used as a pre-conditioner for MCMC change-of-variable proposals.}
  \item{pbd}{Probability of performing a birth/death proposal, otherwise perform a rotate proposal.}
  \item{pb}{Probability of performing a birth proposal given that we choose to perform a birth/death proposal.}
  \item{stepwpert}{Initial width of proposal distribution for peturbing cut-points.}
  \item{probchv}{Probability of performing a change-of-variable proposal.  Otherwise, only do a perturb proposal}.
  \item{minnumbot}{Minimum number of observations required in leaf (terminal) nodes.}
  \item{ntree}{Vector of length two for the number of trees used for the mean
    model and the number of trees used for the variance model.}
  \item{numcut}{Number of cutpoints to use for each predictor variable.}
  \item{xicuts}{More detailed construction of cut-points can be specified
    by the \code{xicuts} function and provided here.}
  \item{power}{Power parameter in the tree depth penalizing prior.}
  \item{base}{Base parameter in the tree depth penalizing prior.}
  \item{k}{Prior hyperparameter for the mean model.}
  \item{sigmaf}{SD of \code{y.train} desired for f function leaf prior.}
\item{dist}{Distribution to be passed to intercept-only AFT model to center \code{y.train}.}
  \item{sigmav}{Initialization of square-root of variance parameter.}
  \item{overalllambda}{A rudimentary estimate of the process standard deviation. Used in calibrating the variance prior.}
  \item{overallnu}{Shape parameter for the variance prior.}
\item{K}{Number of grid points for which to estimate survival probability.}
\item{events}{Grid points for which to estimate survival probability.}
\item{drawMuTau}{Whether to utilize DPM or not.}
  \item{alpha}{Initial value of DPM concentration parameter.}
  \item{alpha.a}{Gamma prior parameter setting for DPM concentration parameter
    where E[\code{alpha}]=\code{alpha.a}/\code{alpha.b}.}
  \item{alpha.b}{See \code{alpha.a} above.}
  \item{alpha.draw}{Whether to draw \code{alpha} or it is fixed at the initial value.}
  \item{neal.m}{The number of additional atoms for Neal 2000 DPM algorithm 8.}
  \item{constrain}{Whether to perform constained DPM or unconstrained.}
  \item{k0.draw}{Whether to fix k0 or draw it if from the DPM LIO prior hierarchy: k0~Gamma(1.5, 7.5).}
  \item{b0.draw}{Whether to fix b0 or draw it from the DPM LIO prior hierarchy: b0~Gamma(2, 1).}
  \item{printevery}{Outputs MCMC algorithm status every printevery iterations.}
}

\details{
  \code{nft()} is the function to fit time-to-event data.  The most general form of the model allowed is
  \eqn{Y({\bf x})=f({\bf x})+s({\bf x})Z}{Y(x)=f(x)+s(x)E}
  where \eqn{E} follows a nonparametric error distribution
  by default.
  %% \eqn{f({\bf x})=\sum_{j=1}^m g({\bf x};T_j,M_j)}{f(x)=sum g(x;T_j,M_j)} and
  %% \eqn{s({\bf x})=\prod_{j=1}^{m^\prime} h({\bf x};T^\prime_j,M^\prime_j)}{s(x)=prod g(x;T'_j,M'_j)},
  %% where the \eqn{g(\cdot;T_j,M_j)}{g(.;T_j,M_j)} represent additive tree components used for modeling the mean and
  %% \eqn{h(\cdot;T^\prime_j,M^\prime_j)}{h(.;T'_j,M'_j)} represent multiplicative tree components used for modeling the variance.

  %% The most common models to fit are a homoscedastic single-tree model, a homoscedastic BART model and a heteroscedastic BART model.

  %% For a BART model, set \code{pbd=c(0.7,0.0)} and \code{ntreeh=1}.  This forces a scalar (homoscedastic) variance term.

  %% For a single-tree model, set \code{pbd=(0.7,0.0)}, \code{ntreeh=1} and \code{ntree=1}. This forces the mean component to be modeled using only one tree.

  %% The heteroscedastic BART model is the default.

}
\value{
  \item{res}{Model fit object of S3 class \code{nft} that is essentially a list.}
}

\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{predict.nft}}
}
\examples{

B=getOption('mc.cores', 1)

data(lung)
str(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1
table(delta)

## this study reports time in days rather than weeks or months
times=lung$time
times=times/7  ## weeks
summary(times)

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])
## lung$sex:        Male=1 Female=2

## token run just to test installation
post=nft(x.train, times, delta, tc=B, K=0,
         nskip=0, ndpost=1, nadapt=1, adaptevery=1)

\dontrun{
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
x.test[ , 2]=rep(1:2, each=N)
K=75
events=seq(0, 150, length.out=K+1)
pred = predict(post, x.test, K=K, events=events[-1],
               XPtr=XPtr, tc=B, FPD=TRUE)

plot(events, c(1, pred$surv.fpd.mean[1:K]), type='l', col=4,
     ylim=0:1, 
     xlab=expression(italic(t)), sub='weeks',
     ylab=expression(italic(S)(italic(t), italic(x))))
lines(events, c(1, pred$surv.fpd.upper[1:K]), lty=2, lwd=2, col=4)
lines(events, c(1, pred$surv.fpd.lower[1:K]), lty=2, lwd=2, col=4)
lines(events, c(1, pred$surv.fpd.mean[K+1:K]), lwd=2, col=2)
lines(events, c(1, pred$surv.fpd.upper[K+1:K]), lty=2, lwd=2, col=2)
lines(events, c(1, pred$surv.fpd.lower[K+1:K]), lty=2, lwd=2, col=2)
legend('topright', c('Adv. lung cancer\nmortality example',
                     'M', 'F'), lwd=2, col=c(0, 4, 2), lty=1)

}
}
