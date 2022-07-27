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
    tc=getOption("mc.cores", 1), 
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
    sigmav=NULL, total.lambda=NA, total.nu=10, mask=NULL,
    ## survival analysis 
    K=100, events=NULL, 
    ## DPM LIO
    drawDPM=1L, 
    alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
    neal.m=2, constrain=1, 
    m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
    a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
    ## misc
    na.rm=FALSE, probs=c(0.025, 0.975), printevery=100
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
  \item{probchv}{Probability of performing a change-of-variable proposal.  Otherwise, only do a perturb proposal.}
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
  \item{total.lambda}{A rudimentary estimate of the process standard deviation. Used in calibrating the variance prior.}
  \item{total.nu}{Shape parameter for the variance prior.}
  \item{mask}{If a proportion is provided, then said quantile
    of \code{max.i s(x.i)} is used to mask non-stationary
  departures (with respect to convergence) above this threshold.}
\item{K}{Number of grid points for which to estimate survival probability.}
\item{events}{Grid points for which to estimate survival probability.}
\item{drawDPM}{Whether to utilize DPM or not.}
  \item{alpha}{Initial value of DPM concentration parameter.}
  \item{alpha.a}{Gamma prior parameter setting for DPM concentration parameter
    where E[\code{alpha}]=\code{alpha.a}/\code{alpha.b}.}
  \item{alpha.b}{See \code{alpha.a} above.}
  \item{alpha.draw}{Whether to draw \code{alpha} or it is fixed at the initial value.}
  \item{neal.m}{The number of additional atoms for Neal 2000 DPM algorithm 8.}
  \item{constrain}{Whether to perform constained DPM or unconstrained.}
  \item{m0}{Center of the error distribution: defaults to zero.}
  \item{k0.a}{First Gamma prior argument for \code{k0}.}
  \item{k0.b}{Second Gamma prior argument for \code{k0}.}
  \item{k0}{Initial value of \code{k0}.}
  \item{k0.draw}{Whether to fix k0 or draw it if from the DPM LIO prior
    hierarchy: \code{k0~Gamma(k0.a, k0.b)}, i.e., \code{E[k0]=k0.a/k0.b}.}
  \item{a0}{First Gamma prior argument for \eqn{tau}.}
  \item{b0.a}{First Gamma prior argument for \code{b0}.}
  \item{b0.b}{Second Gamma prior argument for \code{b0}.}
  \item{b0}{Initial value of \code{b0}.}
  \item{b0.draw}{Whether to fix b0 or draw it from the DPM LIO prior 
    hierarchy: \code{b0~Gamma(b0.a, b0.b)}, i.e.,
    \code{E[b0]=b0.a/b0.b}.}
  \item{na.rm}{Value to be passed to the \code{predict} function.}
  \item{probs}{Value to be passed to the \code{predict} function.}
  \item{printevery}{Outputs MCMC algorithm status every printevery iterations.}
}

\details{
  \code{nft()} is the function to fit time-to-event data.  The most general form of the model allowed is
  \eqn{Y({\bf x})=mu+f({\bf x})+s({\bf x})Z}{Y(x)=mu+f(x)+s(x)E}
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
The \code{nft()} function returns a fit object of S3 class type
  \code{nft} that is essentially a list containing the following items.
}
\value{
  %%\item{res}{Model fit object of S3 class \code{nft} that is
  %%essentially a list.}
  
  \item{ots,oid,ovar,oc,otheta}{These are \code{XPtr}s to the BART
    \eqn{f(x)} objects in RAM that are only available for fits generated
    in the current R session. }

  \item{sts,sid,svar,sc,stheta}{Similarly, these are \code{XPtr}s to the
    HBART \eqn{s(x)} objects.}

  \item{fmu}{The constant \eqn{mu}.}
  
  \item{f.train,s.train}{The trained \eqn{f(x)} and \eqn{s(x)}
  respectively: matrices with \code{ndpost} rows and \eqn{n} columns.}

  \item{f.train.mean,s.train.mean}{The posterior mean of the trained
  \eqn{f(x)} and \eqn{s(x)} respectively: vectors of length \eqn{n}.}

\item{f.trees,s.trees}{Character strings representing the trained fits
  of \eqn{f(x)} and \eqn{s(x)} respectively to facilitate usage of the
  \code{predict} function when \code{XPtr}s are unavailable.}
  
\item{dpalpha}{The draws of the DPM concentration parameter
  \eqn{alpha}.}

\item{dpn,dpn.}{The number of atom clusters per DPM, \eqn{J}, for all
  draws including burn-in and excluding burn-in respectively.}

\item{dpmu}{The draws of the DPM parameter \eqn{mu[i]} where
  \eqn{i=1,...,n} indexes subjects: a matrix with \code{ndpost} rows and
  \eqn{n} columns.}

\item{dpmu.}{The draws of the DPM parameter \eqn{mu[j]} where
  \eqn{j=1,...,J} indexes atom clusters: a matrix with \code{ndpost}
  rows and \eqn{J} columns.}

\item{dpwt.}{The weights for efficient DPM calculations by atom clusters
  (as opposed to subjects) for use with \eqn{dpmu.} (and \eqn{dpsd.};
  see below): a matrix with \code{ndpost} rows and \eqn{J} columns.}
  
\item{dpsd,dpsd.}{Similarly, the draws of the DPM parameter \eqn{tau[i]}
  transformed into the standard deviation \eqn{sigma[i]} for
  convenience.}

\item{dpC}{The indices \eqn{j} for each subject \eqn{i} corresponding to
  their shared atom cluster.}

\item{z.train}{The data values/augmentation draws of \eqn{log t}.}

\item{f.tmind/f.tavgd/f.tmaxd}{The min/average/max tier degree of trees in the \eqn{f} ensemble.}

\item{s.tmind/s.tavgd/s.tmaxd}{The min/average/max tier degree of trees in the \eqn{s} ensemble.}

\item{f.varcount,s.varcount}{Variable importance counts of branch
decision rules for each \eqn{x} of \eqn{f} and \eqn{s} respectively:
matrices with \code{ndpost} rows and \eqn{p} columns.}

  \item{f.varcount.mean,s.varcount.mean}{Similarly, the posterior mean
  of the variable importance counts for each \eqn{x} of \eqn{f} and
  \eqn{s} respectively: vectors of length \eqn{p}.}

  \item{f.varprob,s.varprob}{Similarly, re-weighting the posterior mean
  of the variable importance counts as sum-to-one probabilities for each
  \eqn{x} of \eqn{f} and \eqn{s} respectively: vectors of length
  \eqn{p}.}

\item{LPML}{The log Pseudo-Marginal Likelihood as typically
  calculated for right-/left-censoring.}

\item{pred}{The object returned from the \code{predict} function where
  \code{x.test=x.train} in order to calculate the \code{soffset}
  item that is needed to use \code{predict} when \code{XPtr}s are not
  available. }

\item{soffset}{See \code{pred} above.}
\item{aft}{The AFT model fit used to initialize NFT BART.}
\item{elapsed}{The elapsed time of the run in seconds.}

}

\references{
  Sparapani R., Logan B., Laud P. (2021)
  Nonparametric Failure Time: Time-to-event Machine Learning with Heteroskedastic Bayesian Additive Regression Trees and Low Information Omnibus Dirichlet Process Mixtures
  \emph{MCW Biostatistics Technical Report 72}
  \url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/tr72.pdf?la=en}.
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
         nskip=0, ndpost=10, nadapt=4, adaptevery=1)

\donttest{
file.='lung.rds'
if(file.exists(file.)) {
    post=readRDS(file.)
    XPtr=FALSE
} else {
    set.seed(99)
    post=nft(x.train, times, delta, tc=B, K=0)
    XPtr=TRUE
    ##saveRDS(post, file.)
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
