\name{FPD.wbart}
\title{Friedman's partial dependence function for a previously fitted BART model}
\alias{FPD.wbart}
\alias{FPD.pbart}
\alias{FPD.survbart}
\alias{FPD}
\description{
  Friedman's partial dependence function estimates
  the marginal effect of a variable or variables.
}
\usage{
FPD(object,  ## object returned from BART
    x.test,  ## settings of x.test
    S,       ## indices of subset
    subset.=NULL,
    x.train, 
    probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L,
    seed=99L)
\method{FPD}{wbart}(object, x.test, S,
    subset.=NULL, x.train=object$x.train, probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L, seed = 99L)
\method{FPD}{pbart}(object, x.test, S,
    subset.=NULL, x.train=object$x.train, probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L, seed = 99L)
\method{FPD}{survbart}(object, x.test, S,
    subset.=NULL, x.train=object$tx.test, probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L, seed = 99L)
}
\arguments{
   \item{object}{
     \code{object} returned from previous BART fit.
   }

   \item{x.test}{
     Matrix of Q x L covariate settings (where L is
     the length of \code{S}) in the order of
     columns suggested by \code{S}.
   }

    \item{x.train}{
      Matrix of N x P covariates of the training.
      Typically, automated since it is stored
      as a member of \code{object} by \code{gbart}.
   }

   \item{S}{Indices for the marginal effects subset.}

   \item{subset.}{If \code{S} has one or two dimensions,
     then this function supports the Conditional dependence method,
     i.e., subsetting the data according the grid of values in \code{x.test}.
     For two-dimensions, with a continuous and a category, the
     continuous covariate must be increasing.  If feasible,
     more complex two-dimensional scenarios can be considered in the future.}

   \item{probs}{The lower and upper limits for the credible interval.}
   
   \item{mc.cores}{ Number of threads to utilize. }

   \item{mult.impute}{ For missing data in \code{x.train}, perform
     multiple imputation this many times.}
   \item{seed}{ For missing values in \code{x.train}, perform
     multiple imputation with this seed.}
}
%%% \details{
%%%
%%%}
\value{
  Returns a matrix of predictions corresponding to \code{x.test}.
}
\references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (2001)
   Greedy Function Approximation: a Gradient Boosting Machine
   \emph{The Annals of Statistics}, \bold{29}, 1189--1232.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
}
%%\seealso{
%%  \code{\link{FPD.pbart}}, \code{\link{FPDK.wbart}}
%%}

\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
