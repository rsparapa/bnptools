\name{predict.nft}
\alias{predict.nft}
%\docType{package}
\title{Drawing Posterior Predictive Realizations for NFT BART models.}
\description{
  The function \code{predict.nft()} is the main function for drawing posterior predictive realizations at new inputs using a fitted model stored in a \code{nft} object returned from \code{nft()}.
}
\usage{
\method{predict}{nft}(
            ## data
            object,
            x.test=object$x.train,
            ## multi-threading
            tc=1, ##OpenMP thread count
            ## current process fit vs. previous process fit
            XPtr=TRUE,
            ## predictions
            K=0,
            events=object$events,
            FPD=FALSE,
            probs=c(0.025, 0.975),
            take.logs=TRUE,
            na.rm=FALSE,
            seed=NULL,
            ## default settings for NFT:BART/HBART/DPM
            fmu=object$fmu,
            soffset=object$soffset,
            drawDPM=object$drawDPM,
            ## etc.
            ...)
}
\arguments{
  \item{object}{Object of type \code{nft} from a previous call to \code{nft()}.}
  \item{x.test}{New input settings in the form of a matrix at which to construct predictions.  Defaults to the training inputs.}
  \item{tc}{Number of OpenMP threads to use for parallel computing.}
  \item{XPtr}{If \code{object} was created during the currently running
    R process, then (via an Rcpp \code{XPtr}) the BART/HBART tree
    ensemble objects can be accessed in RAM; otherwise, those objects
    will need to be loaded from their string encodings.}
  \item{K}{The length of the grid of time-points to be used
    for survival predictions.  Set to zero to avoid these calculations
    which can be time-consuming for large data sets.}
  \item{events}{You can specify the grid of time-points; otherwise,
    they are derived from quantiles of the augmented event times.}
  \item{FPD}{Whether to yield the usual predictions or marginal predictions
    calculated by the partial dependence function.}
  \item{probs}{A vector of length two containing the
    lower and upper quantiles to be calculated for the predictions.}
  \item{take.logs}{Whether or not to take logarithms.}
  \item{na.rm}{Whether \code{NA} values should be removed from the
    summaries.}
  \item{seed}{If provided, then this value is used to generate random
    event times from the predictive distribution}.
  \item{fmu}{BART centering parameter for the test data.  Defaults to the value used by \code{nft()} when training the model.}
  \item{soffset}{HBART centering parameter for the test data.  Defaults to the value used by \code{nft()} when training the model.}
   \item{drawDPM}{Whether NFT BART was fit with, or without, DPM.}
   \item{...}{The et cetera objects passed to the \code{predict} method.
     Currently, it has no functionality.}
}
\details{
  \code{predict.nft()} is the main function for calculating posterior predictions and uncertainties once a model has been fit by \code{nft()}.  

  Returns a list with the following entries.
}
\value{
  \item{f.test}{Posterior realizations of the mean function stored in a
    matrix. Omitted if partial dependence functions are performed since
  these will typically be large.}
  \item{s.test}{Posterior realizations of the SD function stored in a matrix. Omitted if partial dependence functions are performed since
  these will typically be large.}
  \item{f.test.mean}{Posterior predictive mean of mean function.}
  \item{f.test.lower}{Posterior predictive lower quantile of mean function.}
  \item{f.test.upper}{Posterior predictive upper quantile of mean function.}
  \item{s.test.mean}{Posterior predictive mean of SD function.}
  \item{s.test.lower}{Posterior predictive lower quantile of SD function.}
  \item{s.test.upper}{Posterior predictive upper quantile of SD
    function.}
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
\code{\link{nft}}
}

