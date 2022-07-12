\name{predict.aftree}
\alias{predict.aftree}
%\docType{package}
\title{Estimating the survival and the hazard for AFT BART models.}
\description{
  The function \code{predict.aftree()} is provided for
  performing posterior inference via test data set estimates
  stored in a \code{aftree} object returned from \code{AFTree()} in a similar
  fashion as that of \code{predict.nft}. N.B.
  the \code{x.test} matrix must be provided on the \code{AFTree()}
  function call.  Here we are only calculating the survival function
  by default, and, if requested, the hazard as well.
}
\usage{
\method{predict}{aftree}(
            ## data
            object,
            ## predictions
            events=NULL,
            FPD=FALSE,
            hazard=FALSE,
            probs=c(0.025, 0.975),
            take.logs=TRUE,
            seed=NULL,
            ## default settings
            ndpost=nrow(object$mix.prop),
            nclust=ncol(object$mix.prop),
            ## etc.
            ...)
}
\arguments{
  \item{object}{Object of type \code{nft} from a previous call to \code{nft()}.}
  \item{events}{You must specify a grid of time-points; however, they
    can be a matrix with rows for each subject.}
  \item{FPD}{Whether to yield the usual predictions or marginal predictions
    calculated by the partial dependence function.}
  \item{hazard}{By default, only the survival function is estimated;
    however, with \code{TRUE}, you can request the hazard as well.}
  \item{probs}{A vector of length two containing the
    lower and upper quantiles to be calculated for the predictions.}
  \item{take.logs}{Whether or not to take logarithms.}
    \item{seed}{If provided, then this value is used to generate random
    natural logarithms of event times from the predictive distribution.}
  \item{ndpost}{The number of MCMC samples generated.}
  \item{nclust}{The number of DPM clusters generated.}
   \item{...}{The et cetera objects passed to the \code{predict} method.
     Currently, it has no functionality.}
}
\details{
  Returns a list with the following entries. If
  \code{hazard=TRUE} is specified, then a similar set of
  entries for the hazard are produced.
}
\value{
  \item{surv.fpd}{Survival function posterior draws on a grid of
    time-points by the partial dependence function when requested.}
  \item{surv.fpd.mean}{Survival function estimates on a grid of
    time-points by the partial dependence function when requested.}
  \item{surv.fpd.lower}{Survival function lower quantiles on a grid of
    time-points by the partial dependence function when requested.}
  \item{surv.fpd.upper}{Survival function upper quantiles on a grid of
    time-points by the partial dependence function when requested.}
}
%% \references{
%%   Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (1998)
%%   Bayesian CART model search. 
%%   \emph{Journal of the American Statistical Association}, \bold{93}, 935--948.

%%   Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (2010)
%%   BART: Bayesian additive regression trees. 
%%   \emph{The Annals of Applied Statistics}, \bold{4}, 266--298.

%%   Pratola, Matthew T. (2016)
%%   Efficient Metropolis Hastings proposal mechanisms for Bayesian regression tree models.
%%   \emph{Bayesian analysis}, \bold{11}, 885--911.
 
%%   Pratola, Matthew T., Chipman, Hugh A., George, Edward I. and McCulloch, Robert E. (2017)
%%   Heteroscedastic BART Using Multiplicative Regression Trees.
%%   \emph{arXiv preprint}, \bold{arXiv:1709.07542}, 1--20.
%% }
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{predict.nft}}
}

