% file MASS/man/Boston.Rd
% copyright (C) 1994-2022 W. N. Venables and B. D. Ripley
% copyright (C) 2024 Rodney Sparapani
%
\name{Boston}
\alias{Boston}
\title{
  Housing Values in Suburbs of Boston: with Corrections and Additions
}
\description{
  The \code{Boston} data frame has 506 rows and 21 columns.
  Original data had 14: new variables are in all-caps.
}
\usage{
data(Boston)
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{OBS.}}{Observation number}
    \item{\code{TOWN}}{Town name}
    \item{\code{TOWN.NO}}{Town number}
    \item{\code{TRACT}}{Census Tract}
    \item{\code{LON}}{Longitude}
    \item{\code{LAT}}{Latitude}
    \item{\code{medv}}{
      median value of owner-occupied homes in $1000s.
    }
    \item{\code{CMEDV}}{
      corrected median value of owner-occupied homes in $1000s.
    }
    \item{\code{crim}}{
      per capita crime rate by town.
    }
    \item{\code{zn}}{
      proportion of residential land zoned for lots over 25,000 sq.ft.
    }
    \item{\code{indus}}{
      proportion of non-retail business acres per town.
    }
    \item{\code{chas}}{
      Charles River dummy variable (= 1 if tract bounds river; 0 otherwise).
    }
    \item{\code{nox}}{
      nitrogen oxides concentration (parts per 10 million).
    }
    \item{\code{rm}}{
      average number of rooms per dwelling.
    }
    \item{\code{age}}{
      proportion of owner-occupied units built prior to 1940.
    }
    \item{\code{dis}}{
      weighted mean of distances to five Boston employment centres.
    }
    \item{\code{rad}}{
      index of accessibility to radial highways.
    }
    \item{\code{tax}}{
      full-value property-tax rate per $10,000.
    }
    \item{\code{ptratio}}{
      pupil-teacher ratio by town.
    }
    \item{\code{black}}{
      \eqn{1000(Bk - 0.63)^2} where \eqn{Bk} is the proportion of blacks
      by town.
    }
    \item{\code{lstat}}{
      lower status of the population (percent).
    }
  }
}
\source{
  Harrison, D. and Rubinfeld, D.L. (1978)
  Hedonic prices and the demand for clean air.
  \emph{J. Environ. Economics and Management}
  \bold{5}, 81--102.

  Belsley D.A., Kuh, E.  and Welsch, R.E. (1980)
  \emph{Regression Diagnostics. Identifying Influential Data and Sources
    of Collinearity.}
  New York: Wiley.

  Gilley, O.W. and R. Kelley Pace (1996)
  On the Harrison and Rubinfeld Data.
  \emph{Journal of Environ. Economics and Management}, \bold{31}, 403--405. 

  Pace, R. Kelley and O.W. Gilley (1997)
  Using the Spatial Configuration of the Data to Improve Estimation.
  \emph{Journal of the Real Estate Finance and Economics}, \bold{14}, 333--340. 
}
\keyword{datasets}
