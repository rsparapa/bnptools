\name{NN.wbart}
\title{Nearest Neighbors marginal for a previously fitted BART model}
\alias{NN.wbart}
\alias{NN.pbart}
\alias{NN.survbart}
\alias{NN}
\description{
  Nearest Neighbors estimates of
  the marginal effect for a variable or variables.
}
\usage{
NN(object,  ## object returned from BART
    x.test,  ## settings of x.test
    S,       ## indices of subset
    nearest,
    x.train, 
    probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L,
    seed=99L)
\method{NN}{wbart}(object, x.test, S,
    nearest, x.train=object$x.train, probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L, seed = 99L)
\method{NN}{pbart}(object, x.test, S,
    nearest, x.train=object$x.train, probs=c(0.025, 0.975),
    mc.cores=getOption('mc.cores', 1L),
    mult.impute=4L, seed = 99L)
\method{NN}{survbart}(object, x.test, S,
    nearest, x.train=object$tx.test, probs=c(0.025, 0.975),
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

   \item{nearest}{The indices of \code{x.train} specifying Nearest Neighbors,
     i.e., subsetting the data according to corresponding the grid of values in
     \code{x.test}.
     }

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

%%\seealso{
%%  \code{\link{NN.pbart}}, \code{\link{NNK.wbart}}
%%}

\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
