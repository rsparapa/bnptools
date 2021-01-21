\name{ucarprice}
\alias{ucarprice}
\docType{data}
\title{
Used Car Prices
}
\description{
Prices of used cars and variables describing the cars.
}
\usage{data("ucarprice")}
\format{
  A data frame with 1,000 observations on the following 7 variables.
  \describe{
    \item{\code{price}}{numeric: price of a car in dollars}
    \item{\code{trim}}{factor: trim of the car,  for example, how the interior is styled}
    \item{\code{isOneOwner}}{factor: has the car had just one owner, t=yes, f=no}
    \item{\code{mileage}}{numeric: number of miles the car has been driven}
    \item{\code{year}}{numeric: model year of the car}
    \item{\code{color}}{factor: color of the car}
    \item{\code{displacement}}{factor: displacement of the car engine}
  }
}
\examples{
data(ucarprice)
}
\keyword{datasets}
