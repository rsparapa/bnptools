\name{predict.nftbart}
\alias{predict.nftbart}
%\docType{package}
\title{Drawing Posterior Predictive Realizations for nftbart models.}
\description{
  The function \code{predict.nftbart()} is the main function for drawing posterior predictive realizations at new inputs using a fitted model stored in a \code{nftbart} object returned from \code{nftbart()}.
}
\usage{
\method{predict}{nftbart}(
object,
x.test=object$x.train,
tc=1,
fmean=mean(object$y.train),
q.lower=0.025,
q.upper=0.975,...)
}
\arguments{
  \item{object}{Object of type \code{nftbart} from a previous call to \code{nftbart()}}
  \item{x.test}{New input settings in the form of an \code{npred x p} matrix at which to construct predictions.  Defaults to the training inputs.}
  \item{tc}{Number of OpenMP threads to use for parallel computing.}
  \item{fmean}{Mean-centering vector for the training data.  Defaults to the value used by \code{nftbart()} when fitting the model.  Usually should be left to the default.}
  \item{q.lower}{Lower quantile to return.}
  \item{q.upper}{Upper quantile to return.}
   \item{...}{not used.}
}
\details{
  \code{predict.nftbart()} is the main function for calculating posterior predictions and uncertainties once a model has been fit by \code{nftbart()}.  

  Returns an object of type \code{nftbart} with the following entries.
}
\value{
  \item{mdraws}{Posterior realizations of the mean function, \eqn{f({\bf x})}{f(x)} stored in an \code{ndpost x npred} matrix.}
  \item{mmean}{Posterior predictive mean of \eqn{f({\bf x})}{f(x)}.}
  \item{smean}{Posterior predictive mean of the standard deviation, \eqn{s({\bf x})}{s(x)}.}
  \item{msd}{Posterior standard deviation of the mean, \eqn{f({\bf x})}{f(x)}.}
  \item{ssd}{Posterior standard deviation of the standard devation, \eqn{s({\bf x})}{s(x)}.}
  \item{m.5}{Posterior median of the mean function realizations, \eqn{f({\bf x})}{f(x)}.}
  \item{m.lower}{Posterior \code{q.lower} quantile of the mean function realizations.}
  \item{m.upper}{Posterior \code{q.upper} quantile of the mean function realizations.}
  \item{s.5}{Posterior median of the standard deviation function realizations, \eqn{s({\bf x})}{s(x)}.}
  \item{s.lower}{Posterior \code{q.lower} quantile of the standard deviation function realizations.}
  \item{s.upper}{Posterior \code{q.upper} quantile of the standard deviation function realizations.}
  \item{q.lower}{Lower quantile used in constructing the above.}
  \item{q.upper}{Upper quantile used in constructing the above.}
}
\references{
  Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (1998)
  Bayesian CART model search. 
  \emph{Journal of the American Statistical Association}, \bold{93}, 935--948.

  Chipman, Hugh A., George, Edward I., and McCulloch, Robert E. (2010)
  BART: Bayesian additive regression trees. 
  \emph{The Annals of Applied Statistics}, \bold{4}, 266--298.

  Pratola, Matthew T. (2016)
  Efficient Metropolis Hastings proposal mechanisms for Bayesian regression tree models.
  \emph{Bayesian analysis}, \bold{11}, 885--911.
 
  Pratola, Matthew T., Chipman, Hugh A., George, Edward I. and McCulloch, Robert E. (2017)
  Heteroscedastic BART Using Multiplicative Regression Trees.
  \emph{arXiv preprint}, \bold{arXiv:1709.07542}, 1--20.
}
\author{
Matthew T. Pratola <mpratola@stat.osu.edu> [aut, cre, cph], Robert E. McCulloch <robert.e.mcculloch@gmail.com> [aut, cre, cph], Hugh Chipman <hugh.chipman@gmail.com> [aut, cph]
Maintainer: Matthew T. Pratola <mpratola@stat.osu.edu>, Robert E. McCulloch <robert.e.mcculloch@gmail.com>
}
\seealso{
\code{\link{nft}}
}

