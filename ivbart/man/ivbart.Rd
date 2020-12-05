\name{ivbart}
\title{IVBART for continuous treatments and outcomes}
\alias{ivbart}

\description{ IVBART is a Bayesian \dQuote{sum-of-trees} model.\cr For a
  numeric treatment \eqn{T}, we have the IV first stage: \eqn{T = f(Z,
  X) + eT} where \eqn{f} follows the BART prior.  For a numeric outcome
  \eqn{Y}, we have the IV second stage: \eqn{Y = \beta T + h(X) + eY}{Y
  = beta T + h(X) + eY} where \eqn{h} follows the BART prior.
  N.B. \eqn{eT} and \eqn{eY} are correlated and their nonparametric
  joint distribution is estimated by Dirichlet Process Mixtures (DPM). }

\usage{
ivbart(
    Z, X, T, Y,
    burn=1000, nd=2000, burnf=1000, burnh=1000,
    keepevery=10,
    betabar=0, Abeta=0.01,
    m=200, nc=100,
    power=2, base=0.95,
    k=2, sigmaf=NA, sigmah=NA,
    v=0.17, nu=2.004, a=0.016,
    n = length(Y),
    Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
    centermeans=TRUE,
    fs = rep(0, n), hs = rep(0, n),
    mTs=0, mYs=0,
    betas=NA, sTs=NA, sYs=NA, gammas=NA,
    printevery=100
)
}
\arguments{
   \item{Z}{ A vector of a single instrument or a matrix of multiple instruments. }
   \item{X}{ A vector of a single confounder or a matrix of multiple
     confounders. If no confounders are observed, then provide a scalar
     constant as the argument such as 0. }
   \item{T}{ A vector of the continuous treatment. }
   \item{Y}{ A vector of the continuous outcome. }
   \item{burn}{ The number of \eqn{\beta}{beta} draws to burn-in and discard. }
   \item{nd}{ The number of \eqn{\beta}{beta} draws to return. }
   \item{burnf}{ The number of \eqn{f} draws to burn-in before the initial value. }
   \item{burnh}{ The number of \eqn{h} draws to burn-in before the initial value. }
   \item{keepevery}{ Every \code{keepevery} draw of \eqn{\beta}{beta} is returned. }
   \item{betabar}{ The mean of \eqn{\beta}{beta} for its Normal prior. }
   \item{Abeta}{ The precision of \eqn{\beta}{beta} for its Normal prior. }
   \item{m}{ The number of trees for the BART prior. }
   \item{nc}{ The number of cut-points for the BART prior. }
   \item{power}{ Power parameter for trees in the BART prior. }
   \item{base}{ Base parameter for tree in the BART prior. }
   \item{k}{ The number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-0.5. }
   \item{sigmaf}{ The SD of \eqn{f}. }
   \item{sigmah}{ The SD of \eqn{h}. }
   \item{v}{ }
   \item{nu}{ }
   \item{a}{ }
   \item{n}{ The number of observations as calculated from \code{length(Y)}
     for convenience. Specifying this parameter on the command line is NOT recommended. }
   \item{Imin}{ The minimum number of DPM clusters. }
   \item{Imax}{ The maximum number of DPM clusters. }
   \item{psi}{ }
   \item{gs}{ }
   \item{centermeans}{ By default, center \code{T} and \code{Y}.  Set this parameter
     to \code{FALSE} to refrain from centering. }
   \item{fs}{ }
   \item{hs}{ }
   \item{mTs}{ The starting value for \eqn{\mu_{Ts}}{muTs}: TSLS is used
     when not provided. }
   \item{mYs}{ The starting value for \eqn{\mu_{Ys}}{muYs}: TSLS is used
     when not provided. }
   \item{betas}{ The starting value for \eqn{\beta}{beta}: TSLS is used
     when not provided. }
   \item{sTs}{ The starting value for \eqn{\sigma_{Ts}}{sdTs}: TSLS is used
     when not provided. }
   \item{sYs}{ The starting value for \eqn{\sigma_{Ys}}{sdYs}: TSLS is used
     when not provided. }
   \item{gammas}{ The starting value for \eqn{\gamma}{gamma}: TSLS is used
     when not provided. }
   \item{printevery}{ During MCMC, a message is printed every \code{printevery} draws. }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case.
   %and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
  \code{ivbart} returns a list containing the following items.
  
   \item{dnpart}{ The number of DPM clusters. }
   \item{dalpha}{ Draws of the DPM concentration parameters. }
   \item{dmu1}{ Draws of the mean from the DPM with a Normal prior for \eqn{eT}. }
   \item{dsigma1}{ Draws of the SD from the DPM with a Normal prior for \eqn{eT}. }
   \item{dmu2}{ Draws of the mean from the DPM with a Normal prior for \eqn{eY}. }
   \item{dsigma2}{ Draws of the SD from the DPM with a Normal prior for \eqn{eY}. }
   \item{dbeta}{ Draws of the \eqn{\beta}{beta} with a Normal prior. }
   \item{df}{ Draws of \eqn{f}. }
   \item{dh}{ Draws of \eqn{h}. }
   \item{dfburn}{ Draws to create an initial value for \eqn{f}. }
   \item{dhburn}{ Draws to create an initial value for \eqn{h}. }
   \item{ag}{ }
   \item{priag}{ }
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.
}
\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu}.
}
\examples{

data(nlsym)

Y=nlsym$lwage76
subset=!is.na(Y)
table(subset)
Y=Y[subset]
T=nlsym$ed76[subset]
Z=cbind(nlsym$nearc2, nlsym$nearc4)[subset, ]
dimnames(Z)[[2]]=c('nearc2', 'nearc4')
X=cbind(nlsym$exp76, (nlsym$exp76-mean(nlsym$exp76))^2,
        nlsym$black, nlsym$smsa76r, nlsym$reg76r)[subset, ]
dimnames(X)[[2]]=c('exp76', 'exp762', 'black',
                   'smsa76r', 'reg76r')

sd.Y = sd(Y)
sd.T = sd(T)
T. = (T-mean(T))/sd.T ## name ending in dot for re-scaled
Y. = (Y-mean(Y))/sd.Y
T = T - mean(T)
Y = Y - mean(Y)
ratio = sd.Y/sd.T

## TSLS
## first stage manually
df.tsls1 = data.frame(T., Z, X)
tsls1 = lm(T.~., df.tsls1)
## second stage manually
df.tsls2 = data.frame(Y., T. = tsls1$fitted.values, X)
tsls2 = lm(Y.~., df.tsls2)

## IVBART
L = t(chol(var(cbind(tsls1$resid, tsls2$resid))))
X. = X[ , -2] ## without the square: BART will not need it

##test IVBART with a token run to ensure installation works
set.seed(99)
post = ivbart(Z, X., T., Y.,
              burn=1, burnf=1, burnh=1, keepevery=1,
              nd=1, sigmaf=1, sigmah=1,
              betas=tsls2$coeff[2], ## beta for scaled data
              sTs=L[1, 1], gammas=L[2, 1], sYs=L[2, 2],
              printevery=1,
              Imin=2, Imax=floor(0.5*length(Y))+1, gs=500)

\dontrun{
    set.seed(99)
    post = ivbart(Z, X., T., Y.,
                  nd=1000, sigmaf=1, sigmah=1,
                  betas=tsls.$coeff[2], ## beta for scaled data
                  sTs=L[1, 1], gammas=L[2, 1], sYs=L[2, 2],
                  Imin=2, Imax=floor(0.5*length(Y))+1, gs=500)

    beta.ivbart = post$dbeta * ratio

    max.=1.2*exp(max(beta.ivbart))
    plot(density(exp(beta.ivbart)),
         xlim=c(1/max., max.), ylim=c(0, 30),
         log='x', lwd=2, main=expression(hat(beta)),
         xlab='Additional year of schooling: multiple of wages')
    abline(h=0)
    abline(v=1, lty=2)
}

}

