\name{tsvs}
\alias{tsvs}
\title{Variable selection with NFT BART models.}
\description{
  The \code{tsvs()} function is for Thompson sampling
  variable selection with NFT BART. 
}
\usage{
tsvs(
               ## data
               x.train, times, delta=NULL, 
               rm.const=TRUE, rm.dupe=TRUE,
               ##tsvs args
               K=20, a.=1, b.=0.5, C=0.5,
               rds.file='tsvs.rds', pdf.file='tsvs.pdf',
               ## multi-threading
               tc=getOption("mc.cores", 1), ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100, 
               chv=NULL,
               method="spearman", use="pairwise.complete.obs",
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(10, 2), numcut=100, xicuts=NULL,
               power=c(2, 2), base=c(0.95, 0.95),
               ## f function
               k=5, sigmaf=NA,
               dist='weibull', 
               ## s function
               sigmav=NULL, total.lambda=NA, total.nu=10, mask=0.95,
               ## survival analysis 
               ##K=100, events=NULL, 
               ## DPM LIO
               drawDPM=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
               a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
               ## misc
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100,
               transposed=FALSE
)

}
\arguments{
  \item{x.train}{nxp matrix of predictor variables for the training data.}
  \item{times}{nx1 vector of the observed times for the training data.}
  \item{delta}{nx1 vector of the time type for the training data:
    0, for right-censoring; 1, for an event; and, 2, for
    left-censoring.}
      \item{rm.const}{To remove constant variables or not.}
    \item{rm.dupe}{To remove duplicate variables or not.}
  \item{K}{The number of Thompson sampling steps to take.
  Not to be confused with the size of the time grid for survival
  distribution estimation.}
  \item{a.}{The prior parameter for successes of a Beta distribution.}
  \item{b.}{The prior parameter for failures of a Beta distribution.}
  \item{C}{The probability cut-off for variable selection.}
  \item{rds.file}{File name to store RDS object containing Thompson
    sampling parameters.}
  \item{pdf.file}{File name to store PDF graphic of variables selected.}
%\item{x.test}{mxp matrix of predictor variables for the test set.}
%  \item{impute.bin}{Indices of the columns of \code{x.train} to be imputed.}
%  \item{impute.prob}{nx1 vector of prior probabilities for imputation.}
  \item{tc}{Number of OpenMP threads to use.}
  \item{nskip}{Number of MCMC iterations to burn-in and discard.}
  \item{ndpost}{Number of MCMC iterations kept after burn-in.}
  \item{nadapt}{Number of MCMC iterations for adaptation prior to burn-in.}
  \item{adaptevery}{Adapt MCMC proposal distributions every
    \code{adaptevery} iteration.}
    \item{chv}{Predictor correlation matrix used as a pre-conditioner for MCMC change-of-variable proposals.}
    \item{method,use}{Correlation options for change-of-variable proposal
    pre-conditioner.}
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
%\item{K}{Number of grid points for which to estimate survival probability.}
%\item{events}{Grid points for which to estimate survival probability.}
%\item{TSVS}{Setting to \code{TRUE} will avoid unnecessary processing for
%  Thompson sampling variable selection, i.e., all that is needed is
%  the variable counts from the tree branch decision rules.}
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
  \item{printevery}{Outputs MCMC algorithm status every printevery
    iterations.}
      \item{transposed}{\code{tsvs} handles all of the pre-processing
	for \code{x.train/x.test} (including
    tranposing) computational efficiency.}
}

\details{
  \code{tsvs()} is the function to perform variable selection.
  %% \eqn{f({\bf x})=\sum_{j=1}^m g({\bf x};T_j,M_j)}{f(x)=sum g(x;T_j,M_j)} and
  %% \eqn{s({\bf x})=\prod_{j=1}^{m^\prime} h({\bf x};T^\prime_j,M^\prime_j)}{s(x)=prod g(x;T'_j,M'_j)},
  %% where the \eqn{g(\cdot;T_j,M_j)}{g(.;T_j,M_j)} represent additive tree components used for modeling the mean and
  %% \eqn{h(\cdot;T^\prime_j,M^\prime_j)}{h(.;T'_j,M'_j)} represent multiplicative tree components used for modeling the variance.
  %% The most common models to fit are a homoscedastic single-tree model, a homoscedastic BART model and a heteroscedastic BART model.
  %% For a BART model, set \code{pbd=c(0.7,0.0)} and \code{ntreeh=1}.  This forces a scalar (homoscedastic) variance term.
  %% For a single-tree model, set \code{pbd=(0.7,0.0)}, \code{ntreeh=1} and \code{ntree=1}. This forces the mean component to be modeled using only one tree.
  %% The heteroscedastic BART model is the default.
The \code{tsvs()} function returns a fit object of S3 class type
\code{list} as well as storing it in \code{rds.file} for
sampling in progress.
}


\references{
  Sparapani R., Logan B., Laud P. (2021)
  Nonparametric Failure Time: Time-to-event Machine Learning with Heteroskedastic Bayesian Additive Regression Trees and Low Information Omnibus Dirichlet Process Mixtures
  \emph{MCW Biostatistics Technical Report 72}
  \url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/tr72.pdf?la=en}.

  Liu Y., Rockova V. (2021)
  Variable selection via Thompson sampling.
  \emph{Journal of the American Statistical Association. Jun 29:1-8.}
}

\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{tsvs2}}
}
\examples{

##library(nftbart)
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

\donttest{
##vars=tsvs(x.train, times, delta)
vars=tsvs(x.train, times, delta, K=0) ## K=0 just returns 0
}
}
