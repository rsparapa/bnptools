\name{SHAP.wbart}
\title{Shapley values for a previously fitted BART model}
\alias{SHAP.wbart}
\alias{SHAP}
\alias{SHAP2.wbart}
\alias{SHAP2}
\alias{SHNN.wbart}
\alias{SHNN}
\alias{SHNN2.wbart}
\alias{SHNN2}
\description{
  Shapley values are additive estimates for
  the marginal effect of a variable or variables.
  \code{SHAP} for simple random sampling
  and \code{SHNN} for nearest neighbor sampling.
  \code{SHAP2} and \code{SHNN2} for 2-way interactions.
}
\usage{
SHAP(object,  ## object returned from BART
     x.test,  ## settings of x.test
     S,       ## indices of subset
     x.train,
     probs,
     seed,
     mult.impute)
\method{SHAP}{wbart}(object, x.test, S, x.train=object$x.train,
              probs=c(0.025, 0.975), seed = 99L, mult.impute=5L)
SHNN(object,  ## object returned from BART
     x.test,  ## settings of x.test
     S,       ## indices of subset
     x.train,
     probs,
     seed,
     mult.impute)
\method{SHNN}{wbart}(object, x.test, S, x.train=object$x.train,
              probs=c(0.025, 0.975), seed = 99L, mult.impute=5L)
SHAP2(object, ## object returned from BART
     x.test,  ## settings of x.test
     S,       ## indices of subset
     x.train,
     probs,
     seed,
     mult.impute)
\method{SHAP2}{wbart}(object, x.test, S, x.train=object$x.train,
              probs=c(0.025, 0.975), seed = 99L, mult.impute=5L)
SHNN2(object, ## object returned from BART
     x.test,  ## settings of x.test
     S,       ## indices of subset
     x.train,
     probs,
     seed,
     mult.impute)
\method{SHNN2}{wbart}(object, x.test, S, x.train=object$x.train,
              probs=c(0.025, 0.975), seed = 99L, mult.impute=5L)
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

   \item{probs}{The lower and upper limits for the credible interval.}
   
   \item{mult.impute}{Perform
     multiple imputation this many times: default is \code{5}.}
   \item{seed}{Perform
     multiple imputation with this seed: default is \code{99}.}
}
%%% \details{
%%%
%%%}
\value{
  Returns a matrix of predictions corresponding to \code{x.test}.
}
%\references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

%% Friedman, J.H. (2001)
%%    Greedy Function Approximation: a Gradient Boosting Machine
%%    \emph{The Annals of Statistics}, \bold{29}, 1189--1232.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%}
%%\seealso{
%%  \code{\link{SHAP.pbart}}, \code{\link{SHAPK.wbart}}
%%}

\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
