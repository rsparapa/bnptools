\name{ptrtbart}
\title{Additive Treatment BART for dichotomous outcomes with Normal latents}
\alias{ptrtbart}
\description{
  BART is a Bayesian \dQuote{sum-of-trees} model.\cr
  %% For numeric response \eqn{y}, we have
  %% \eqn{y = f(x) + \epsilon}{y = f(x) + e},
  %% where \eqn{\epsilon \sim N(0, 1)}{e ~ N(0, 1)}.\cr
  For a binary response \eqn{y}, the additive treatment BART model is \eqn{P(Y=1 | x) = F(\beta A + f(x))}, where \eqn{F}
  denotes the standard Normal CDF (probit link) and $A$ denotes the treatment indicator.
  
  In this case, \eqn{f} is the sum of many tree models.
  The goal is to have very flexible inference for the unknown
  function \eqn{f} while adjusting for the treatment using a parametric model.
  
  In the spirit of \dQuote{ensemble models}, each tree is constrained by a prior to be a weak learner
  so that it contributes a  small amount to the overall fit.
}
\usage{
  ptrtbart(
  x.train1,trt.train, y.train, 
  x.test1=matrix(0.0,0,0),trt.test=0,
  mu0=0,sig20=25,
  sparse=FALSE, theta=0, omega=1,
  a=0.5, b=1, augment=FALSE, rho=NULL,
  x1info=matrix(0.0,0,0), usequants=FALSE,
  cont=FALSE, rm.const=TRUE,
  k=2.0, power=2.0, base=.95,
  binaryOffset=NULL,
  ntree=200L, numcut1=100L, 
  ndpost=1000L, nskip=100L, keepevery=1L,
  nkeeptrain=ndpost, nkeeptest=ndpost,
  ##nkeeptestmean=ndpost,
  nkeeptreedraws=ndpost,
  printevery=100L, transposed=FALSE
  ##treesaslists=FALSE
)
}
\arguments{
  \item{x.train1}{
    Explanatory variables for training (in sample) data without the treatment.\cr
    May be a matrix or a data frame,
    with (as usual) rows corresponding to observations and columns to variables.\cr
    If a variable is a factor in a data frame, it is replaced with dummies.
    Note that q dummies are created if q>2 and
    one dummy is created if q=2, where q is the number of levels of the factor.
    \code{ptrtbart} will generate draws of \eqn{f(x)} for each \eqn{x}
    which is a row of \code{x.train1}.
  }
  \item{trt.train}{ Treatment variable for training (in sample) data. \cr It is given as a vector.
  }
  \item{y.train}{
    Binary dependent variable for training (in sample) data.\cr
    %% If y is numeric a continous response model is fit (Normal errors).\cr
    %% If y is a factor (or just has values 0 and 1) then a binary response model
    %% with a probit link is fit.
  }
  
  \item{x.test1}{
    Explanatory variables for test (out of sample) data without the treatment.\cr
    Should have same structure as x.train1.\cr
    \code{ptrtbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test1}.
  }
  \item{trt.test}{Treatment variable for testing(out of sample) data. Should have same structure as trt.train. }
  \item{mu0}{ Prior mean for the treatment effect. Default value is zero}
 \item{sig20}{ Prior variance for the treatment effect. Default value is 25.}
  \item{sparse}{Whether to perform variable selection based on a
    sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
  \item{theta}{Set \eqn{theta} parameter; zero means random.}
  \item{omega}{Set \eqn{omega} parameter; zero means random.}
  \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
      \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
  \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
    \eqn{b=1}.}
  \item{rho}{Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the
    number of covariates under consideration.}
  \item{augment}{Whether data augmentation is to be performed in sparse
    variable selection.}
  
  \item{x1info}{ You can provide the cutpoints to BART or let BART
    choose them for you.  To provide them, use the \code{x1info}
    argument to specify a list (matrix) where the items (rows) are the
    covariates and the contents of the items (columns) are the
    cutpoints.  }
  
  \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }
  
  \item{cont}{ Whether or not to assume all variables are continuous.}
  
  \item{rm.const}{ Whether or not to remove constant variables.}
  
  \item{k}{
    %% For numeric y,
    %% k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
    %% The response (y.train) is internally scaled to range from -.5 to .5.
    For binary y,
    k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
    %In both cases,
    The bigger k is, the more conservative the fitting will be.
  }
  
  \item{power}{
    Power parameter for tree prior.
  }
  
  \item{base}{
    Base parameter for tree prior.
  }
  
  \item{binaryOffset}{
    Used for binary \eqn{y}.\cr
    The model is \eqn{P(Y=1 | x) = F( \beta A + f(x) + binaryOffset)}.\cr
    %% The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
    %% a probability other than .5.
  }
  
  %% \item{M}{
    %%   For large data sets, to efficiently sample \eqn{f} with Modified LISA MCMC
    %%   divide your data into \eqn{M} shards.
    %%   Defaults to 1 shard, i.e., plain BART.  See Entezari et al.
    %% }
  
  \item{ntree}{
    The number of trees in the sum.
  }
  
  \item{numcut1}{
    The number of possible values of c (see usequants).
    If a single number if given, this is used for all variables.
    Otherwise a vector with length equal to ncol(x.train1) is required,
    where the \eqn{i^{th}}{i^th} element gives the number of c used for
    the \eqn{i^{th}}{i^th} variable in x.train1.
    If usequants is false, numcut equally spaced cutoffs
    are used covering the range of values in the corresponding
    column of x.train1.  If usequants is true, then  min(numcut, the number of unique values in the
                                                        corresponding columns of x.train1 - 1) c values are used.
  }
  
  \item{ndpost}{
    The number of posterior draws returned.
  }
  
  \item{nskip}{
    Number of MCMC iterations to be treated as burn in.
  }
  
  \item{nkeeptrain}{
    Number of MCMC iterations to be returned for train data.
  }
  
  \item{nkeeptest}{
    Number of MCMC iterations to be returned for test data.
  }
  
  %% \item{nkeeptestmean}{
    %% Number of MCMC iterations to be returned for test mean.
    %% }
  
  \item{nkeeptreedraws}{
    Number of MCMC iterations to be returned for tree draws.
  }
  
  \item{keepevery}{
    Every keepevery draw is kept to be returned to the user.
  }
  
  \item{printevery}{
    As the MCMC runs, a message is printed every printevery draws.
  }
  
  %% \item{keepevery}{
    %% Every keepevery draw is kept to be returned to the user.\cr
    %% A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
    %% and \eqn{f^*(x)}{f*(x)}
    %% at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
    %% the current draw of \eqn{f}.
    %% }
  
  \item{transposed}{
    When running \code{ptrtbart} in parallel, it is more memory-efficient
    to transpose \code{x.train1} and \code{x.test1}, if any, prior to
    calling \code{mc.ptrtbart}.
  }
  
  %% \item{treesaslists}{
    %%   The default behavior is to return the trees as a text string.  If you
    %%   specify \code{treesaslists=TRUE}, then the trees will be returned
    %%   as a list as well.
    %% }
}
\details{
  BART is an Bayesian MCMC method.
  At each MCMC interation, we produce a draw from %% the joint posterior
  %% \eqn{f | (x,y)}{f \| (x,y)} in the numeric \eqn{y} case
  %% and just
  \eqn{f} in the binary \eqn{y} case.
  
  Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
  from which fits and summaries may be extracted.  The output consists of values
  \eqn{f^*(x)}{f*(x)}
  %(and \eqn{\sigma^*}{sigma*} in the numeric case)
  where * denotes a particular draw.
  The \eqn{x} is either a row from the training data (x.train1) or the test data (x.test).
}
\value{
  \code{ptrtbart} returns an object of type \code{ptrtbart} which is
  essentially a list. % assigned class \sQuote{bart}.
  
  \item{yhat.train1}{ A matrix with ndpost rows and nrow(x.train1) columns.
    Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{ f}
    and each column corresponds to a row of x.train1.
    The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
    and the \eqn{j^{th}}{j\^th} row of x.train1.\cr
    Burn-in is dropped.
  }
  
  \item{beta}{ Post burn in draws of beta, the treatment effect, length=ndpost}
  
 \item{yhat.train1}{ A matrix with ndpost rows and nrow(x.train1) columns.
    Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{ f}
    and each column corresponds to a row of x.train1.
    The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
    and the \eqn{j^{th}}{j\^th} row of x.train1.\cr
    Burn-in is dropped.
  }
  
  \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
    %   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
    %   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
    \item{varcount}{a matrix with ndpost rows and nrow(x.train1) columns.
    Each row is for a draw. For each variable (corresponding to the columns),
    the total count of the number of times
    that variable is used in a tree decision rule (over all trees) is given.}
    
    In addition the list
    has a binaryOffset component giving the value used.
    
    Note that in the binary \eqn{y}, case yhat.train and yhat.test are
    \eqn{ \beta A + f(x) } + binaryOffset.  If you want draws of the probability
    \eqn{P(Y=1 | x)} you need to apply the Normal CDF (\code{pnorm})
    to these values.
  }
    \references{
    Chipman, H., George, E., and McCulloch R. (2010)
    Bayesian Additive Regression Trees.
    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi: 10.1214/09-AOAS285>.
    
    Chipman, H., George, E., and McCulloch R. (2006)
    Bayesian Ensemble Learning.
    Advances in Neural Information Processing Systems 19,
    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
    
    Entezari, R., Craiu, R.V., and Rosenthal, J.S. (2017)
    Likelihood Inflating Sampling Algorithm.
    \url{https://arxiv.org/pdf/1605.02113.pdf}.
    
    Friedman, J.H. (1991)
    Multivariate adaptive regression splines.
    \emph{The Annals of Statistics}, \bold{19}, 1--67.
    
    Linero, A.R. (2016)
    Bayesian regression trees for high dimensional prediction and variable
    selection. \emph{JASA},
    \url{http://dx.doi.org/10.1080/01621459.2016.1264957}  
    }
    \author{
    Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
    Rodney Sparapani: \email{rsparapa@mcw.edu}.
    }
    \seealso{
    \code{\link{wbart}}
    }
    \examples{
    
    #Generate data from a probit model
    library(BART)
      N <- 1000
      gX <- NULL; fZ <- NULL; hX <- NULL
      X1 <- rnorm(N)
      X2 <- rnorm(N)
      X3 <- rnorm(N) 
      A <- rbinom(N,1,0.5)
      fZ <- A
      gX <- 3 + 3*X1 - 3*X2 - 2*X3
      hX <- -2*X1 - 7*X2 - X3 + 3*X2*X3 
      
      #Center the above functions
      
      gX <- gX - mean(gX)
      hX <- hX - mean(hX)
      fZ <- fZ - mean(fZ)
      lY <- 0.2368*gX + 2.0125*fZ + 0*fZ*hX  #The interaction effect if zero since the model is additive
      pr <- pnorm(lY)
      U  <-  runif(N,0,1)
      Y  <-  ifelse(U < pr,1,0)  
      
      data <- data.frame(TRT=A,X1=X1,X2=X2,X3=X3,Y=Y,pr=pr) 
      xtrain1 <-  data[, grep('X',colnames(data))]
      TRT <-  data[, grep('TRT',colnames(data))]  
      ytrain <- data$Y
      
      # Fit an additive treatment BART model 
      ndpost <- 1000
      post_ptrtbart <- ptrtbart(x.train1 = xtrain1,
                             trt.train =  TRT,
                             y.train = ytrain,
                             ndpost = ndpost) 
                             
     plot(density(post_ptrtbart$beta),xlab=expression(beta),
     main="Posterior distribution for the treatment effect (beta)")
     abline(v=2.0125,col="red") #Red line representes the true treatment effect 
     
     #Obtain posterior samples for the fit from BART
     
     fitbart <-  post_ptrtbart$yhat.train1
     
     #Obtain posterior samples for the fit from the parametric model for the treatment
     
     fittrt <- matrix(NA,nrow=ndpost,ncol=length(TRT))
     for(i in 1:length(TRT)){
               for(d in 1:ndpost){
                  fittrt[d,i] <-  post_ptrtbart$beta[d]*TRT[i]   
   }
   }
   
   #Add the fitted values from each component and the centering parameter; 
   #Average across the posterior samples
   
   fitadd <-  colMeans(fitbart + fittrt + post_ptrtbart$binaryOffset)
   
   #The combined fitted values is also stored as yhat.train
   
   fitall <- colMeans(post_ptrtbart$yhat.train)
   
   #verify that the two methods give the same results
   
  plot(fitall,fitadd) 
   
   #We can also verity that the fitted probabilities are the same by applying the standard normal CDF (pnorm) to the posterior samples of the fitted values.
   
   
    }
    }
    \keyword{nonparametric}
    \keyword{tree}
    \keyword{additive models}
    \keyword{regression}
    \keyword{nonlinear}
    