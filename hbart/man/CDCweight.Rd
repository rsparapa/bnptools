\name{CDCweight}
\alias{CDCweight}
\docType{data}
\title{CDC weight for age growth charts}

\description{ Using the Cole and Green LMS method, here we provide
  percentiles of weight by age and sex based on the US National Center for Health
  Statistics data for children aged 2 to 17. }

\usage{
data(CDCweight)
}

\format{
  \tabular{ll}{
    age:\tab Age in years\cr
    sex:\tab 1=male, 2=female\cr
    weight.XXX:\tab Weight XXXth percentile in kg\cr
    M:\tab mean\cr
    S:\tab coefficient of variation\cr
    L:\tab Box-Cox power transformation\cr
    SD:\tab SD=mean X coefficient of variation\cr
  }
}

\references{
  Cole, Timothy J and Green, Pamela J (1992)
  Smoothing reference centile curves: the LMS method and penalized
  likelihood. \emph{Statistics in medicine}, \bold{11}, 1305--1319.

  The US Centers for Disease Control and Prevention stature by age LMS parameters
  \url{https://www.cdc.gov/growthcharts/data/zscore/wtage.csv}
}

\keyword{datasets}
