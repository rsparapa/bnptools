\name{polychor}
\alias{polychor}

\title{Polychoric Correlation}
\description{
  Computes the polychoric correlation (and its standard error)
  between two ordinal variables or from their contingency table, under the 
  assumption that the ordinal variables dissect continuous latent variables that are bivariate normal. Either
  the maximum-likelihood estimator or a (possibly much) quicker ``two-step'' approximation is available. For the ML
  estimator, the estimates of the thresholds and the covariance matrix of the estimates are also available.
}

\usage{
polychor(x, y, ML = FALSE, control = list(), 
  std.err = FALSE, maxcor=0.9999, start, thresholds=FALSE)
}

\arguments{
  \item{x}{a contingency table of counts or an ordered categorical variable; the latter can be numeric, logical, a factor,
    an ordered factor, or a character variable, but if a factor, its levels 
    should be in proper order, and the values of a character variable are
    ordered alphabetically.}
  \item{y}{if \code{x} is a variable, a second ordered categorical variable.}
  \item{ML}{if \code{TRUE}, compute the maximum-likelihood estimate; if \code{FALSE}, the default, compute a quicker
    ``two-step'' approximation.}
  \item{control}{optional arguments to be passed to the \code{optim} function.}
  \item{std.err}{if \code{TRUE}, return the estimated variance of the correlation (for the two-step estimator) 
    or the estimated covariance matrix (for the ML estimator) of the correlation and thresholds; the default is \code{FALSE}.}
  \item{maxcor}{maximum absolute correlation (to insure numerical stability).}
    \item{start}{optional start value(s): if a single number, start value for the correlation; if a list with the elements \code{rho}, \code{row.thresholds}, and \code{column.thresholds}, start values for these parameters; start values are supplied automatically if omitted, and are only relevant when the ML estimator or standard errors are selected.}
    \item{thresholds}{if \code{TRUE} (the default is \code{FALSE}) return estimated thresholds along with the estimated correlation even if standard errors aren't computed.}
}

\value{
  If \code{std.err} or \code{thresholds} is \code{TRUE}, 
  returns a list %an object of class \code{"polycor"}
  with the following components:
  \item{type}{set to \code{"polychoric"}.}
  \item{rho}{the polychoric correlation.}
  \item{row.cuts}{estimated thresholds for the row variable (\code{x}), for the ML estimate.}
  \item{col.cuts}{estimated thresholds for the column variable (\code{y}), for the ML estimate.}
  \item{var}{the estimated variance of the correlation, or, for the ML estimate, 
    the estimated covariance matrix of the correlation and thresholds.}
  \item{n}{the number of observations on which the correlation is based.}
  \item{chisq}{chi-square test for bivariate normality.}
  \item{df}{degrees of freedom for the test of bivariate normality.}
  \item{ML}{\code{TRUE} for the ML estimate, \code{FALSE} for the two-step estimate.}
  Othewise, returns the polychoric correlation.
}

\details{
  The ML estimator is computed by maximizing the bivariate-normal likelihood with respect to the 
  thresholds for the two variables (\eqn{\tau^{x}_i, i = 1,\ldots, r - 1}{\tau^x[i], i = 1,\ldots, r - 1}; 
  \eqn{\tau^{y}_j, j = 1,\ldots, c - 1}{\tau^y[j], j = 1,\ldots, c - 1}) and
  the population correlation (\eqn{\rho}). Here, \eqn{r} and \eqn{c} are respectively the number of levels
  of \eqn{x} and \eqn{y}. The likelihood is maximized numerically using the \code{\link{optim}} function,
  and the covariance matrix of the estimated parameters is based on the numerical Hessian computed by \code{optim}.
  
  The two-step estimator is computed by first estimating the thresholds (\eqn{\tau^{x}_i, i = 1,\ldots, r - 1}{\tau^x[i], i = 1,\ldots, r - 1} 
  and \eqn{\tau^{y}_j, i = j,\ldots, c - 1}{\tau^y[j], i = j,\ldots, c - 1}) separately from the marginal distribution of each variable. Then the 
  one-dimensional likelihood for \eqn{\rho} is maximized numerically, using \code{\link{optim}} if standard errors are
  requested, or \code{\link{optimise}} if they are not. The standard error computed treats the thresholds as fixed.
}

\references{
    Drasgow, F. (1986) 
    Polychoric and polyserial correlations. 
    Pp. 68--74 in S. Kotz and N. Johnson, eds.,
    \emph{The Encyclopedia of Statistics, Volume 7.} Wiley.
    
    Olsson, U. (1979)
    Maximum likelihood estimation of the polychoric correlation coefficient.
    \emph{Psychometrika} \bold{44}, 443-460.
}

\author{John Fox}

\seealso{\code{\link{polyserial}}, \code{\link[stats]{optim}}}

\examples{
if(require(mvtnorm)){
    set.seed(12345)
    data <- rmvnorm(1000, c(0, 0), matrix(c(1, .5, .5, 1), 2, 2))
    x <- data[,1]
    y <- data[,2]
    cor(x, y)  # sample correlation
    }
if(require(mvtnorm)){
    x <- cut(x, c(-Inf, .75, Inf))
    y <- cut(y, c(-Inf, -1, .5, 1.5, Inf))
    polychor(x, y)  # 2-step estimate
    }
if(require(mvtnorm)){
    polychor(x, y, ML=TRUE, std.err=TRUE)  # ML estimate
    }
}

%\keyword{models}
