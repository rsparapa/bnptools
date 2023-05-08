\name{crisk.pre.bart}
\title{Data construction for competing risks with BART}
\alias{crisk.pre.bart}
\description{
Competing risks contained in \eqn{(t, \delta, x)} must be translated to data
suitable for the BART competing risks model; see \code{crisk.bart} for more details.
}
\usage{
crisk.pre.bart( times, delta, x.train=NULL, x.test=NULL,
                x.train2=x.train, x.test2=x.test, K=NULL,
                rm.const=TRUE, numcut=100, grp=NULL,
                xinfo=matrix(0, 0, 0), usequants=FALSE )
}
\arguments{
   \item{times}{
     The time of event or right-censoring.\cr
   }

   \item{delta}{
     The event indicator: 1 is a cause 1 event, 2 a cause 2 while 0 is censored.\cr
   }

   \item{x.train}{
   Explanatory variables for training (in sample) data of cause 1.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data of cause 1.\cr
   If provided, must be a matrix and have the same structure as x.train.\cr
 }
 
   \item{x.train2}{
   Explanatory variables for training (in sample) data of cause 2.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test2}{
   Explanatory variables for test (out of sample) data of cause 2.\cr
   If provided, must be a matrix and have the same structure as x.train.\cr
 }
 
   \item{K}{
     If provided, then coarsen \code{times} per the quantiles
     \eqn{1/K, 2/K, ..., K/K}.
   }

      \item{xinfo}{ You can provide the cutpoints to BART or let BART
     choose them for you.  To provide them, use the \code{xinfo}
     argument to specify a list (matrix) where the items (rows) are the
     covariates and the contents of the items (columns) are the
     cutpoints.  }

      \item{numcut}{ The number of possible values of \eqn{c} (see
   \code{usequants}).  If a single number if given, this is used for all
   variables.  Otherwise a vector with length equal to
   \code{ncol(x.train)} is required, where the \eqn{i^{th}}{i^th}
   element gives the number of \eqn{c} used for the \eqn{i^{th}}{i^th}
   variable in \code{x.train}.  If usequants is false, numcut equally
   spaced cutoffs are used covering the range of values in the
   corresponding column of \code{x.train}.  If \code{usequants} is true, then
   \eqn{min(numcut, the number of unique values in the corresponding
     columns of x.train - 1)} values are used.  }
 
   \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }
   
  \item{rm.const}{ Whether or not to remove constant variables.}

  \item{grp}{A vector of inverse weights for the Dirichlet prior arguments.
     If all the variables are continuous, then \code{grp} is a vector of 1s.
     However, for categorical variables (like factors in a data.frame), the
   inverse weights are the number of categories.  See \code{bartModelMatrix}
   for the details of the default automated derivation when \code{grp=NULL}.
 }
 
}

\value{
   \code{surv.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses for cause 1.}
   \item{y.train2}{A vector of binary responses for cause 2.}
   \item{cond}{A vector of indices of \code{y.train} indicating censored subjects.}
   \item{binaryOffset}{The binary offset for \code{y.train}.}
   \item{binaryOffset2}{The binary offset for \code{y.train2}.}
   \item{tx.train}{A matrix with rows consisting of time and the
   covariates of the training data for cause 1.}
   \item{tx.train2}{A matrix with rows consisting of time and the
   covariates of the training data for cause 2.}
   \item{tx.test}{A matrix with rows consisting of time and the
   covariates of the test data, if any, for cause 1.}
   \item{tx.test2}{A matrix with rows consisting of time and the
   covariates of the test data, if any, for cause 2.}
}

%% \references{
%%   Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%    \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
%% }
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{crisk.bart}}
}
\examples{

data(transplant)

delta <- (as.numeric(transplant$event)-1)

delta[delta==1] <- 4
delta[delta==2] <- 1
delta[delta>1] <- 2
table(delta, transplant$event)

table(1+floor(transplant$futime/30.5)) ## months
times <- 1+floor(transplant$futime/30.5)

x.train=data.frame(abo=transplant[ , 3])
i=which(transplant$abo=='O')[1]
i[2]=which(transplant$abo=='A')[1]
i[3]=which(transplant$abo=='B')[1]
i[4]=which(transplant$abo=='AB')[1]
x.test=data.frame(abo=transplant[i, 3])

pre <- crisk.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

}
\keyword{data construction}
