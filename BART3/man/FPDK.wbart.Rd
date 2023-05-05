\name{FPDK.wbart}
\title{Friedman's partial dependence function with kernel sampling for a previously fitted BART model}
\alias{FPDK.wbart}
%\alias{predict.pbart}
%\alias{predict.lbart}
%\alias{predict.mbart}
\description{
  Friedman's partial dependence function estimates
   with kernel sampling for the marginal effect of a variable or variables.
}
\usage{
\method{FPDK}{wbart}(object, x.test,
    S, x.train=object$x.train,
    mult.impute=4L, kern.var=TRUE,
    alpha=0.05, probs=c(0.025, 0.975), 
    mc.cores=getOption('mc.cores', 1L),
    seed = 99L, nice=19L)
}
\arguments{
   \item{object}{
     \code{object} returned from previous BART fit.
   }

   \item{x.test}{
     Matrix of Q x P covariate settings. Only the
     subset \code{x.test[ , S]} are actually
     used.  Nevertheless, all P columns must be provided.
   }

    \item{x.train}{
      Matrix of N x P covariates of the training.
      Typically, automated since it is stored
      as a member of \code{object} by \code{gbart}.
   }

   \item{S}{Indices for the marginal effects subset.}

   \item{mc.cores}{ Number of threads to utilize. }
   \item{kern.var}{Whether to perform the kernel sampling variance adjustment.}
   \item{alpha}{ Kernel sampling symmetric credible interval level. }
   \item{probs}{ Kernel sampling credible interval limits. }
   \item{mult.impute}{ The number of kernel samples to perform.}
   \item{seed}{ For kernel sampling of \code{x.train}, perform
     multiple imputation with this seed.}
   \item{nice}{ Niceness priority for threads on \code{"unix"}.}
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

Lundberg, S.M. and Lee, S.I. (2017)
   A unified approach to interpreting model predictions
   \emph{Advances in neural information processing systems}, \bold{30}, 1--10.

Janzing, D. and Minorics, L. and Blobaum, P. (2020)
   Feature relevance quantification in explainable AI: A causal problem
   \emph{International conference on artificial intelligence and
     statistics}, 2907--2916.
   
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
