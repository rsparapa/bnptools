\name{gbart}
\title{Generalized BART for continuous and binary outcomes}
\alias{gbart}
\alias{mc.gbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For a numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma^2)}.\cr

\eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a small amount to the overall fit.
}
\usage{
gbart(
      x.train, y.train,
      x.test=matrix(0,0,0), type='wbart',
      ntype=as.integer(
          factor(type, levels=c('wbart', 'pbart', 'lbart'))),
      sparse=FALSE, theta=0, omega=1,
      a=0.5, b=1, augment=FALSE, rho=0, grp=NULL, varprob=NULL,
      xinfo=matrix(0,0,0), usequants=FALSE,
      rm.const=TRUE,
      sigest=NA, sigdf=3, sigquant=0.90,
      k=2, power=2, base=0.95,
      impute.mult=NULL, impute.prob=NULL,
      impute.miss=NULL,
      %sigmaf=NA,
      lambda=NA, tau.num=c(NA, 3, 6)[ntype], %tau.interval=0.9973,
      offset=NULL, w=rep(1, length(y.train)),
      ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
      %ntree=200L, numcut=100L,
      ndpost=1000L, nskip=100L, %keepevery=1L,
      keepevery=c(1L, 10L, 10L)[ntype],
      printevery=100L, transposed=FALSE,
      probs=c(0.025, 0.975),
      mc.cores = 1L, ## mc.gbart only
      nice = 19L,    ## mc.gbart only
      seed = 99L,    ## mc.gbart only
      meta = FALSE,  ## mc.gbart only
      TSVS = FALSE,  ## gbart only
      verbose = 1L, shards = 1L, weight=rep(NA, shards)
)

mc.gbart(
         x.train=matrix(0,0,0), y.train=NULL,
         x.test=matrix(0,0,0), type='wbart',
         ntype=as.integer(
             factor(type, levels=c('wbart', 'pbart', 'lbart'))),
         sparse=FALSE, theta=0, omega=1,
         a=0.5, b=1, augment=FALSE, rho=0, grp=NULL, varprob=NULL,
         xinfo=matrix(0,0,0), usequants=FALSE,
         rm.const=TRUE,
         sigest=NA, sigdf=3, sigquant=0.90,
         k=2, power=2, base=0.95,
         impute.mult=NULL, impute.prob=NULL,
         impute.miss=NULL,
         %sigmaf=NA,
         lambda=NA, tau.num=c(NA, 3, 6)[ntype], %tau.interval=0.9973,
         offset=NULL, w=rep(1, length(y.train)),
         %ntree=200L, numcut=100L,
         ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
         ndpost=1000L, nskip=100L, %keepevery=1L,
         keepevery=c(1L, 10L, 10L)[ntype],
         printevery=100L, transposed=FALSE,
         probs=c(0.025, 0.975),
         mc.cores = getOption('mc.cores', 2L), ## mc.gbart only
         nice = 19L,  ## mc.gbart only
         seed = 99L,  ## mc.gbart only
         meta = FALSE,## mc.gbart only
         TSVS = FALSE,  ## gbart only
         verbose = 1L, shards = 1L, weight=rep(NA, shards)
)

}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
    data.\cr May be a matrix or a data frame, with (as usual) rows
    corresponding to observations and columns to variables.\cr If a
    variable is a factor in a data frame, it is replaced with dummies.
    Note that \eqn{q} dummies are created if \eqn{q>2} and one dummy
    created if \eqn{q=2} where \eqn{q} is the number of levels of the
    factor.  \code{gbart} will generate draws of \eqn{f(x)} for each
    \eqn{x} which is a row of \code{x.train}.  }

   \item{y.train}{
   Continuous or binary dependent variable for training (in sample) data.\cr
If \eqn{y} is numeric, then a continuous BART model is fit (Normal errors).\cr
If \eqn{y} is binary (has only 0's and 1's), then a binary BART model
with a probit link is fit by default: you can over-ride the default via the
argument \code{type} to specify a logit BART model.
   }

   \item{x.test}{ Explanatory variables for test (out of sample)
   data. Should have same structure as \code{x.train}.
   \code{gbart} will generate draws of \eqn{f(x)} for each \eqn{x} which
   is a row of \code{x.test}.  }

 \item{type}{ You can use this argument to specify the type of fit.
   \code{'wbart'} for continuous BART, \code{'pbart'} for probit BART or
   \code{'lbart'} for logit BART. }

 \item{ntype}{ The integer equivalent of \code{type} where
  \code{'wbart'} is 1, \code{'pbart'} is 2 and
  \code{'lbart'} is 3.}
 %\item{rfinit}{ Whether to initialize BART with a greedy RandomForest
 %  fit: the default is \code{FALSE}.}   
   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
   \item{theta}{Set \eqn{theta} parameter; zero means random.}
   \item{omega}{Set \eqn{omega} parameter; zero means random.}
   \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
     \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
   \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
     \eqn{b=1}.}

   \item{augment}{Whether data augmentation is to be performed in sparse
     variable selection.}
   \item{rho}{A multiplier for the inverse weights of the Dirichlet prior
     arguments.  For sparsity, \eqn{rho=p} where \eqn{p} is the
     number of covariates under consideration: the default, \code{rho=0}
     means \eqn{rho=p} (computed by \code{rho=sum(1/grp)}). For more sparsity,
     \eqn{rho<p}, set this argument manually.  See also \code{grp}.}
   \item{grp}{A vector of inverse weights for the Dirichlet prior arguments.
     If all the variables are continuous, then \code{grp} is a vector of 1s.
     However, for categorical variables (like factors in a data.frame), the
   inverse weights are the number of categories.  See \code{bartModelMatrix}
   for the details of the default automated derivation when \code{grp=NULL}.
     }   
   \item{varprob}{ You initialize the variable selection probability:
     defaults to \code{NULL} that means \eqn{1/p}.}
   
   \item{xinfo}{ You can provide the cutpoints to BART or let BART
     choose them for you.  To provide them, use the \code{xinfo}
     argument to specify a list (matrix) where the items (rows) are the
     covariates and the contents of the items (columns) are the
     cutpoints.  }

   \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }
   
   \item{rm.const}{ Whether or not to remove constant variables.}
  
   \item{sigest}{ The prior for the error variance
   (\eqn{sigma^2}{sigma\^2}) is inverted chi-squared (the standard
   conditionally conjugate prior).  The prior is specified by choosing
   the degrees of freedom, a rough estimate of the corresponding
   standard deviation and a quantile to put this rough estimate at.  If
   \code{sigest=NA} then the rough estimate will be the usual least squares
   estimator.  Otherwise the supplied value will be used.
   Not used if \eqn{y} is binary.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
   Not used if \eqn{y} is binary.
   }

   \item{sigquant}{ The quantile of the prior that the rough estimate
   (see \code{sigest}) is placed at.  The closer the quantile is to 1, the more
   aggresive the fit will be as you are putting more prior weight on
   error standard deviations (\eqn{sigma}) less than the rough
   estimate.  Not used if \eqn{y} is binary.  }

   \item{k}{ For numeric \eqn{y}, \code{k} is the number of prior
     standard deviations \eqn{E(Y|x) = f(x)} is away from +/-0.5.
     %The response, code{y.train}, is internally scaled to range from
     %-0.5 to 0.5.  
     For binary \eqn{y}, \code{k} is the number of prior standard
   deviations \eqn{f(x)} is away from +/-3.  The bigger \code{k} is, the more
   conservative the fitting will be.  }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
 }
 \item{impute.mult}{ A vector of the columns of \code{x.train}
   which are multinomial indicators that require imputation:
   the default is \code{NULL}.}
 \item{impute.prob}{ A matrix of probabilities for the
   multinomial indicators that require imputation:
   the default is \code{NULL}.}
 \item{impute.miss}{ A vector of missing indicators for
   the multinomial indicators that require imputation:
   the default is \code{NULL}.}
   %% \item{sigmaf}{
   %% The SD of \eqn{f}.  Not used if \eqn{y} is binary.
   %% }

   \item{lambda}{
   The scale of the prior for the variance.  If \code{lambda} is zero,
     then the variance is to be considered fixed and known at the given
     value of \code{sigest}.  Not used if \eqn{y} is binary.
 }
 
 \item{tau.num}{ The numerator in the \code{tau} definition, i.e.,
   \code{tau=tau.num/(k*sqrt(ntree))}. }
   %% \item{tau.interval}{
   %%   The width of the interval to scale the variance for the terminal
   %%   leaf values.  Only used if \eqn{y} is binary.}

   \item{offset}{ Continous BART operates on \code{y.train} centered by
   \code{offset} which defaults to \code{mean(y.train)}.  With binary
   BART, the centering is \eqn{P(Y=1 | x) = F(f(x) + offset)} where
   \code{offset} defaults to \code{F^{-1}(mean(y.train))}.  You can use
   the \code{offset} parameter to over-ride these defaults.}
   
   \item{w}{ Vector of weights which multiply the standard deviation.
   Not used if \eqn{y} is binary.  }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{numcut}{ The number of possible values of \eqn{c} (see
   \code{usequants}).  If a single number if given, this is used for all
   variables.  Otherwise a vector with length equal to
   \code{ncol(x.train)} is required, where the \eqn{i^{th}}{i^th}
   element gives the number of \eqn{c} used for the \eqn{i^{th}}{i^th}
   variable in \code{x.train}.  If usequants is false, numcut equally
   spaced cutoffs are used covering the range of values in the
   corresponding column of \code{x.train}.  If \code{usequants} is true, then
   \eqn{min(numcut, the number of unique values in the corresponding
   columns of x.train - 1)} values are used.  }

   \item{ndpost}{
   The number of posterior draws returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   %% A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   %% and \eqn{f^*(x)}{f*(x)}
   %% at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   %% the current draw of \eqn{f}.
   }

   \item{transposed}{
   When running \code{gbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.gbart}.}
   \item{probs}{ The lower and upper quantiles to summarize:
     the default is \code{c(0.025, 0.975)}.}

 %% \item{hostname}{
 %%   When running on a cluster occasionally it is useful
 %%   to track on which node each chain is running; to do so
 %%   set this argument to \code{TRUE}.
 %% }
 
    \item{seed}{
     Setting the seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     Number of cores to employ in parallel.
   }

   \item{nice}{
     Set the job niceness.  The default
     niceness is 19: niceness goes from 0 (highest) to 19 (lowest).
   }
   \item{TSVS}{ If \code{TRUE}, then avoid unnecessary calculations
     for speed: \code{gbart} only.}
   \item{verbose}{ If set to \code{0L}, then compute silently.}
   \item{shards}{ For the Modified LISA method, this is the number of
     shards: the default is \code{1L}.}
   \item{weight}{ For the Modified LISA method, this is a vector of
     weights to combine the shards: the default is \code{rep(NA, shards)}.}
   \item{meta}{ Whether or not to produce meta-analysis-like
     estimates from a sharded analysis (as opposed to a Modified LISA
     approach): default is \code{FALSE}.}
}
\details{
   BART is a Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce
   a single model object from which fits and summaries may be extracted.
   The output consists of values \eqn{f^*(x)}{f*(x)} (and
   \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a
   particular draw.  The \eqn{x} is either a row from the training data,
   \code{x.train} or the test data, \code{x.test}.

   For \code{x.train}/\code{x.test} with missing data elements, \code{gbart}
   will singly impute them with hot decking. For one or more missing
   covariates, record-level hot-decking imputation \cite{deWaPann11} is
   employed that is biased towards the null, i.e., nonmissing values
   from another record are randomly selected regardless of the
   outcome. Since \code{mc.gbart} runs multiple \code{gbart} threads in
   parallel, \code{mc.gbart} performs multiple imputation with hot
   decking, i.e., a separate imputation for each thread.  This
   record-level hot-decking imputation is biased towards the null, i.e.,
   nonmissing values from another record are randomly selected
   regardless of \code{y.train}.

 }

 \value{
   %% The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   %% The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   %% including the burn-in draws.  Initially these draws will decline as BART finds fit
   %% and then level off when the MCMC has burnt in.\cr
   %% The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   %% the corresponding \eqn{f(x)} on the vertical axis.

   \code{gbart} returns an object of type \code{gbart} which is
   essentially a list. % assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

   \item{offset}{The data centering value for the BART prior.}
   \item{x.train}{The training data returned with any updates
     due to missing value imputation, factor expansion, etc.}
   \item{yhat.train}{
   A matrix with \code{ndpost} rows and \code{nrow(x.train)} columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }
   \item{yhat.test}{Same as yhat.train but now with \code{x.test} data.}
   \item{yhat.*.mean}{mean of \code{yhat.train/test} fit.}
   \item{yhat.*.lower}{lower quantile of \code{yhat.train/test} fit.}
   \item{yhat.*.upper}{upper quantile of \code{yhat.train/test} fit.}
   \item{prob.train/test}{produced for dichotomous outcomes.}
   \item{prob.*.mean}{mean of \code{prob.train/test} fit.}
   \item{prob.*.lower}{lower quantile of \code{prob.train/test} fit.}
   \item{prob.*.upper}{upper quantile of \code{prob.train/test} fit.}
   \item{sigma}{all draws of sigma including burn-in.}
   \item{sigma.}{only kept draws of sigma with burn-in discarded.}
   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
 }
 \item{treedraws}{A list containing the tree draws.}
   \item{varcount}{matrix with \code{ndpost} rows and \code{nrow(x.train)} columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}
   \item{varprob}{instead of counts, this is the probability that each
   variable is chosen.}
 \item{var*.mean}{The \code{varcount/prob} mean of its posterior.}
 \item{accept}{the accept probability from Metropolis-Hastings step
   within Gibbs.}
 \item{chains}{The number of MCMC chains.}
 \item{grp}{The grouping variable for grouped variables with the
   DART sparse prior.}
 \item{proc.time}{The time elapsed as returned by \code{proc.time()}.}
 \item{LPML}{The Log Pseudo-Marginal Likelihood.  Beware for
   nonparametric models like BART, this quantity can be unstable.}
}
%% \references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

%% De Waal, T., Pannekoek, J. and Scholtus, S. (2011)
%%    Handbook of statistical data editing and imputation.
%%    John Wiley & Sons, Hoboken, NJ.
  
%% Friedman, J.H. (1991)
%%    Multivariate adaptive regression splines.
%%    \emph{The Annals of Statistics}, \bold{19}, 1--67.

%% Gramacy, RB and Polson, NG (2012)
%%    Simulation-based regularized logistic regression.
%%    \emph{Bayesian Analysis}, \bold{7}, 567--590.

%% Holmes, C and Held, L (2006)
%%    Bayesian auxiliary variable models for binary and multinomial regression.
%%    \emph{Bayesian Analysis}, \bold{1}, 145--68.

%% Linero, A.R. (2016)
%%   Bayesian regression trees for high dimensional prediction and variable
%%   selection. \emph{JASA},
%%   \url{http://dx.doi.org/10.1080/01621459.2016.1264957}  
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{bartModelMatrix}}
}
\examples{
##simulate data (example from Friedman MARS paper)
f = function(x){
10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
n = 100      ##number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) ##10 variables, only first 5 matter
y=f(x)+rnorm(n)

##test BART with token run to ensure installation works
set.seed(99)
bartFit = gbart(x,y,nskip=5,ndpost=5)

\dontrun{
##run BART
set.seed(99)
bartFit = gbart(x,y)
}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
