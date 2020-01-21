\name{gbmm}
\title{Generalized BART Mixed Model for continuous and binary outcomes}
\alias{gbmm}
\alias{mc.gbmm}
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

This mixed BART model uses a Bayesian parametric random
effect term \eqn{za} where \eqn{a \sim N(0,\Sigma_a)}
and \eqn{z} is the random effect design vector. This term is combined
with the Bayesian non-parametric BART model for
the fixed effects as written above. This function currently estimates
continuous outcomes \eqn{y} only.
}
\usage{
gbmm(
      x.train, y.train,
      x.test=matrix(0,0,0), type='wbart',
      id.train=NULL, z.train=NULL,
      ntype=as.integer(
          factor(type, levels=c('wbart', 'pbart'))),
      sparse=FALSE, theta=0, thetaDraw=1, omega=1,
      a=0.5, b=1, augment=FALSE, rho=NULL,
      xinfo=matrix(0,0,0), usequants=FALSE,
      rm.const=TRUE,
      sigest=NA, sigdf=3, sigquant=0.90,
      k=2, power=2, base=0.95,
      lambda=NA, tau.num=c(NA, 3, 6)[ntype],
      mxps=NULL,
      offset=NULL, 
      ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
      ndpost=1000L, nskip=100L, %keepevery=1L,
      keepevery=c(1L, 10L, 10L)[ntype],
      printevery=100L, transposed=FALSE,
      hostname=FALSE,
      mc.cores = 1L, ## mc.gbmm only
      nice = 19L,    ## mc.gbmm only
      seed = 99L     ## mc.gbmm only
)

mc.gbmm(
         x.train, y.train,
         x.test=matrix(0,0,0), type='wbart',
         id.train=NULL, z.train=NULL,
         ntype=as.integer(
             factor(type, levels=c('wbart', 'pbart'))),
         sparse=FALSE, theta=0, thetaDraw=1, omega=1,
         a=0.5, b=1, augment=FALSE, rho=NULL,
         xinfo=matrix(0,0,0), usequants=FALSE,
         rm.const=TRUE,
         sigest=NA, sigdf=3, sigquant=0.90,
         k=2, power=2, base=0.95,
         lambda=NA, tau.num=c(NA, 3, 6)[ntype],
         mxps=NULL,
         offset=NULL, 
         ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
         ndpost=1000L, nskip=100L, %keepevery=1L,
         keepevery=c(1L, 10L, 10L)[ntype],
         printevery=100L, transposed=FALSE,
         hostname=FALSE,
         mc.cores = 2L, nice = 19L, seed = 99L
)

}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
    data.\cr May be a matrix or a data frame, with (as usual) rows
    corresponding to observations and columns to variables.\cr If a
    variable is a factor in a data frame, it is replaced with dummies.
    Note that \eqn{q} dummies are created if \eqn{q>2} and one dummy
    created if \eqn{q=2} where \eqn{q} is the number of levels of the
    factor.  \code{gbmm} will generate draws of \eqn{f(x)} for each
    \eqn{x} which is a row of \code{x.train}.  }

   \item{y.train}{
   Continuous variable for training (in sample) data.\cr
If \eqn{y} is numeric, then a continuous BART model is fit (Normal errors).\cr
   }

   \item{x.test}{ Explanatory variables for test (out of sample)
   data. Should have same structure as \code{x.train}.
   \code{gbmm} will generate draws of \eqn{f(x)} for each \eqn{x} which
   is a row of \code{x.test}.  }
 
   \item{id.train}{
   For a two level model, \code{id.train} should be a vector of group
   indices. For a three or more level model, \code{id.train} should be
   either a list of vectors (one for each level) of the group indices
   for each level or a matrix with each column of the matrix specifying
   the group indices for each level.
   }

   \item{z.train}{
     Generalized random effect design matrix. The default is a
     vector of repeating \code{1}'s, indicating a random intercept
     model. A single matrix is required for the two level model, while a
     list of matrices is required for the three or more level
     model. Each matrix should have as many rows as y.train has
     elements. Note that a singleton \code{1} passed as an element of
     the list (or \code{z.train} itself for a two level model) is a shortcut
     for a random intercept model.     
   }
   
   \item{type}{ You can use this argument to specify the type of outcome.
     \code{'wbart'} for continuous outcomes and \code{'pbart'} for
     probit link binary outcomes.
   Note: currently \code{'cgbmm'} only allows continuous outcomes so this
   parameter has no effect.}

 \item{ntype}{ The integer equivalent of \code{type} where
  \code{'wbart'} is 1 and \code{'pbart'} is 2. Note: currently \code{'cgbmm'}
  only allows continuous outcomes so this parameter has no effects.}

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
   \item{theta}{Set \eqn{theta} concentration parameter for Dirichlet;
     zero means random. }
   \item{thetaDraw}{How to sample theta (1=gridding;2=slice sampling)}
   \item{omega}{Set \eqn{omega} parameter; zero means random. Currently
     not used.}
   \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
     \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
   \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
     \eqn{b=1}.}
   \item{rho}{Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the
     number of covariates under consideration.}
   \item{augment}{Whether data augmentation is to be performed in sparse
     variable selection.}
   
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
   standard deviations \eqn{E(Y|x) = f(x)} is away from +/-0.5.  The
   response, code{y.train}, is internally scaled to range from -0.5 to
   0.5.  For binary \eqn{y}, \code{k} is the number of prior standard
   deviations \eqn{f(x)} is away from +/-3.  The bigger \code{k} is, the more
   conservative the fitting will be.  }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

 \item{mxps}{
   A list of prior settings for the mixed effect variance
   components. For the two level model, there should be three elements, scalar
   \code{prior} which specifes the prior distribution (1=Hier.
   Inverse-Wishart (default); 2=Inverse-Wishart), scalar \code{df} which specifies the degree of
   freedom parameter (default of 3), and vector \code{scale} which specifes the scale
   parameters (default of 1). For a three or more level model, this should be specified
   as a list, one element for each level, where each element has the
   same components as described above.
   }

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
   When running \code{gbmm} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.gbmm}.
 }

 \item{hostname}{
   When running on a cluster occasionally it is useful
   to track on which node each chain is running; to do so
   set this argument to \code{TRUE}.
 }
 
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

   For \code{x.train}/{x.test} with missing data elements, \code{gbmm}
   will singly impute them with hot decking. For one or more missing
   covariates, record-level hot-decking imputation \cite{deWaPann11} is
   employed that is biased towards the null, i.e., nonmissing values
   from another record are randomly selected regardless of the
   outcome. Since \code{mc.gbmm} runs multiple \code{gbmm} threads in
   parallel, \code{mc.gbmm} performs multiple imputation with hot
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

   \code{gbmm} returns an object of type \code{gbmm} which is
   essentially a list. % assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.train}{
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}
   \item{varprob}{a matrix with ndpost rows and nrow(x.train)
   columns. Each row is for a draw. For each variable (corresponding to
   the columns), the Dirichlet draw of the variability splitting
   probability is given.
   }
   \item{re.train}{a matrix with the rows representing the posterior
     draws and the columns representing the different groups.
   }
   \item{re.varcov}{an array representing the posterior draws of the
   variance-covariance matrix, where the first dimension represents the
   posterior draws and the second and third dimension represent the
   variance-covariance matrix.}
   \item{re.train.mean}{re.train averaged across the posterior
     distribution.}
   \item{re.varcov.mean}{re.varcov averaged across the posterior
     distribtion.}
   \item{re.corr}{re.varcov transformed to a correlation matrix.}
   \item{re.corr.mean}{re.corr averaged across the posterior distribution.}
   \item{sigest}{
   The rough estimate of the standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

De Waal, T., Pannekoek, J. and Scholtus, S. (2011)
   Handbook of statistical data editing and imputation.
   John Wiley & Sons, Hoboken, NJ.
  
Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
   \emph{The Annals of Statistics}, \bold{19}, 1--67.

Linero, A.R. (2018)
  Bayesian regression trees for high dimensional prediction and variable selection.
  \emph{JASA}, \bold{113(522)}, 626--636.
}
\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu},\cr
Charles Spanbauer: \email{cspanbauer@mcw.edu}.
}
\examples{
library(mxBART)
## Covariates
x <- matrix(runif(10000),ncol=10)
## Random effects
u <- rep(rnorm(100),each=10)
v <- rep(rnorm(50),each=20)
## Second level labels
id2 <- rep(1:100,each=10)
## Third level labels
id3 <- rep(1:50,each=20)
## Outcome; u+v means the random intercept 3 level model is correct
y <- sin(pi*x[,1]*x[,2])+2*(x[,3]-.5)^2+x[,4]+.5*x[,5]+u+v+rnorm(1000,sd=.5)

## Random intercept two level model (no z.train needs to be passed)
ri2.fit <- gbmm(y.train=y,
                x.train=x,
                id.train=id2)

## Random slope two level model
rs2.fit <- gbmm(y.train=y,
                x.train=x,
                id.train=id2,
		z.train=cbind(1,x[,1]))

## Random intercept three level model; CORRECT MODEL!!! 
ri3.fit <- gbmm(y.train=y,
                x.train=x,
                id.train=cbind(id2,id3))
mean(ri3.fit$sigma)
ri3.fit$re.varcov.mean

## Same as previous, but with z.train specified explicitly
ri3b.fit <- gbmm(y.train=y,
                x.train=x,
                id.train=cbind(id2,id3),
                z.train=list(1,1))

## Random slope three level model
## (random intercept on level two and random slope on third level)
rs3.fit <- gbmm(y.train=y,
                x.train=x,
                id.train=cbind(id2,id3),
                z.train=list(1,cbind(1,x[,1])))

## Passing different settings
## Flat prior (large scale parameters as opposed to default of 1)
flat.fit <- gbmm(y.train=y,
                 x.train=x,
                 id.train=cbind(id2,id3),
                 mxps=list(list(prior=1,df=1,scale=1000),
                           list(prior=1,df=1,scale=1000)))
}

\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{random effects}
\keyword{nonlinear}
