\name{galaxy}
\alias{galaxy}
\docType{data}
\title{Galaxy velocities}
\description{
  This data set consider physical information on velocities (km/second) for
  82 galaxies reported by Roeder (1990). These are drawn from six well-separated
  conic sections of the Corona Borealis region.
}
\usage{data(galaxy)}
\format{
  A data frame with 82 observations on the following variable.
  \describe{
    \item{\code{speed}}{a numeric vector giving the speed of galaxies ((km/second))}
  }
}
\source{
  Roeder, K. (1990) Density estimation with confidence sets exemplified by 
     superclusters and voids in the galaxies, Journal of the American Statistical
     Association, 85: 617-624.
}
\references{

  Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
     Using Mixtures. Journal of the American Statistical Association, 90: 577-588.
  
}
\examples{
data(galaxy)
## maybe str(galaxy) ; plot(galaxy) ...
}
\keyword{datasets}
