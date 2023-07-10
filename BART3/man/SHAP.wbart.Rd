\name{SHAP.wbart}
\title{Shapley values for a previously fitted BART model}
\alias{SHAP.wbart}
%\alias{predict.pbart}
%\alias{predict.lbart}
%\alias{predict.mbart}
\description{
  Shapley values are additive estimates for
  the marginal effect of a variable or variables.
}
\usage{
\method{SHAP}{wbart}(object, x.test, S, x.train=object$x.train,
              type='wbart', probs=c(0.025, 0.975), call=FALSE)
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
   
   \item{call}{ Whether to use C++ or R code: FALSE for R code, TRUE for
     C+.}

   \item{type}{The type of BART object fit.}

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
