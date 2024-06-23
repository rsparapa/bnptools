\name{maxRhat}
\alias{maxRhat}
\alias{oldRhat}
\alias{splitRhat}
\title{Rhat for MCMC convergence testing}
\description{
  Rhat has evolved over several generations. \code{maxRhat()}
  is the latest and greatest.  We also provide \code{oldRhat()}
  for comparison.  And the most recent prior generation is an
  important component: \code{splitRhat()}.
}

\usage{
maxRhat(y, C)
splitRhat(y, C)
oldRhat(y, C)
}

\arguments{
\item{y}{A vector of parameters drawn from \code{C}
    MCMC posterior chains combined by simple concatenation
    where burn-in was discarded.
}
\item{C}{The number of chains.}
}

\details{ The cutoff for \code{oldRhat} was considered to be 1.1 while
  for \code{maxRhat} it was lowered to 1.01.  However, these diagnostics
  were developed for parametric models.  In our experience with BART,
  \code{maxRhat} values so low are practically non-existent.
  Unfortunately, there does
  not appear to be any research in this area for Bayesian nonparametric
  models.  So we have adopted a rule of thumb for the cutoff, 1.1, that
  appears to be adequate for convergence checking with BART.  Of course,
  more work needs to be done particularly with respect to BART.  }

\value{
  A list that includes at least \code{Rhat}.
  \item{minSeff}{An estimate of the Effective Sample Size provided by
  \code{maxRhat} only.}
  \item{minSpct}{An estimate of the Effective Sample Size as a
    percentage provided by \code{maxRhat} only.}
}

\references{
  A Vehtari, A Gelman, D Simpson, B Carpenter and PC Burkner (2021).
  Rank-normalization, folding and localization: an improved Rhat for
  assessing convergence of MCMC.
  \emph{Bayesian Analysis}, \bold{16}, 667-718.
}

\seealso{
   \code{\link{gewekediag}}
}
 
\keyword{convergence diagnostics}

\examples{

\dontrun{
##simulate data (example from Friedman MARS paper)
f = function(x){
10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
n = 100      ##number of observations
set.seed(21)
x=matrix(runif(n*10),n,10) ##10 variables, only first 5 matter
y=f(x)+rnorm(n)

##run BART
C <- 8
options(mc.cores = C)
bartFit = mc.gbart(x,y,seed = 99)

print(maxRhat(bartFit$sigma., C))
}
}
