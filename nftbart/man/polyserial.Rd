\name{polyserial}
\alias{polyserial}

\title{Polyserial Correlation}

\description{
  Computes the polyserial correlation (and its standard error)
  between a quantitative variable and an ordinal variable, based on the assumption
  that the joint distribution of the quantitative variable and a latent continuous variable underlying the ordinal 
  variable is bivariate normal. Either
  the maximum-likelihood estimator or a quicker ``two-step'' approximation is available. For the ML
  estimator the estimates of the thresholds and the covariance matrix of the estimates are also available.
}
\usage{
polyserial(x, y, ML = FALSE, control = list(), 
  std.err = FALSE, maxcor=.9999, bins=4, start, thresholds=FALSE)
}

\arguments{
  \item{x}{a numerical variable.}
  \item{y}{an ordered categorical variable; can be numeric, logical, a factor,
   an ordered factor, or a character variables,
   but if a factor, its levels should be in proper order, and the values of a
   character variable are ordered alphabetically.}
  \item{ML}{if \code{TRUE}, compute the maximum-likelihood estimate; if \code{FALSE}, the default, compute a quicker
    ``two-step'' approximation.}
  \item{control}{optional arguments to be passed to the \code{optim} function.}
  \item{std.err}{if \code{TRUE}, return the estimated variance of the correlation (for the two-step estimator)
    or the estimated covariance matrix of the correlation and thresholds (for the ML estimator);
    the default is \code{FALSE}.}
  \item{maxcor}{maximum absolute correlation (to insure numerical stability).}
  \item{bins}{the number of bins into which to dissect \code{x} for a test of bivariate normality; the default is 4.}
    \item{start}{optional start value(s): if a single number, start value for the correlation; if a list with the elements \code{rho} and \code{thresholds}, start values for these parameters; start values are supplied automatically if omitted, and are only relevant when the ML estimator or standard errors are selected.}
    \item{thresholds}{if \code{TRUE} (the default is \code{FALSE}) return estimated thresholds along with the estimated correlation even if standard errors aren't computed.}
}

\value{
  If \code{std.err} or \code{thresholds} is \code{TRUE}, 
  returns a list %an object of class \code{"polycor"}
  with the following components:
  \item{type}{set to \code{"polyserial"}.}
  \item{rho}{the polyserial correlation.}
  \item{cuts}{estimated thresholds for the ordinal variable (\code{y}), for the ML estimator.}
  \item{var}{the estimated variance of the correlation, or, for the ML estimator, 
    the estimated covariance matrix of the correlation and thresholds.}
  \item{n}{the number of observations on which the correlation is based.}
  \item{chisq}{chi-square test for bivariate normality.}
  \item{df}{degrees of freedom for the test of bivariate normality.}
  \item{ML}{\code{TRUE} for the ML estimate, \code{FALSE} for the two-step estimate.}
  Othewise, returns the polyserial correlation.
}

\details{
  The ML estimator is computed by maximizing the bivariate-normal likelihood with respect to the 
  thresholds for \eqn{y} (\eqn{\tau^{y}_j, i = 1,\ldots, c - 1}{\tau^y[j], j = 1,\ldots, c - 1}) and
  the population correlation (\eqn{\rho}). The likelihood is maximized numerically using the \code{\link{optim}} function,
  and the covariance matrix of the estimated parameters is based on the numerical Hessian computed by \code{optim}.
  
  The two-step estimator is computed by first estimating the thresholds 
  (\eqn{\tau^{y}_j, i = 1,\ldots, c - 1}{\tau^y[j], j = 1,\ldots, c - 1})
  from the marginal distribution of \eqn{y}. Then if the standard error of \eqn{\hat{\rho}}{\rho hat} is requested, the 
  one-dimensional likelihood for \eqn{\rho} is maximized numerically, using \code{\link{optim}} if standard errors are
  requested; the standard error computed treats the thresholds as fixed. If the standard error isn't requested,
  \eqn{\hat{\rho}}{\rho hat} is computed directly. 
}

\references{
    Drasgow, F. (1986) 
    Polychoric and polyserial correlations. 
    Pp. 68--74 in S. Kotz and N. Johnson, eds.,
    \emph{The Encyclopedia of Statistics, Volume 7.} Wiley.
}

\author{John Fox}

\seealso{\code{\link{polychor}}, \code{\link[stats]{optim}}}

\examples{
if(require(mvtnorm)){
    set.seed(12345)
    data <- rmvnorm(1000, c(0, 0), matrix(c(1, .5, .5, 1), 2, 2))
    x <- data[,1]
    y <- data[,2]
    cor(x, y)  # sample correlation
    }
if(require(mvtnorm)){
    y <- cut(y, c(-Inf, -1, .5, 1.5, Inf))
    polyserial(x, y)  # 2-step estimate
    }
if(require(mvtnorm)){
    polyserial(x, y, ML=TRUE, std.err=TRUE) # ML estimate
    }
}

%\keyword{models}

