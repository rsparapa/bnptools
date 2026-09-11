
\name{polyserial}
\alias{polyserial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for calculating polyserial correlations  }
\description{
 This function calculates polyserial correlations using either an ad hoc or ML estimator.
}
\usage{
polyserial(x, y, ml = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A continuous variable. }
  \item{y}{ An ordinal variable with at least two categories. }
  \item{ml}{A logical variable indicating whether to use a formal maximum likelihood estimator (default) or an ad hoc estimator. }
}
\details{
The variables should be numeric. The function returns NA if y has only one category.}
\value{ Returns the polyserial correlation. }
\references{ Olsson, U., Drasgow, F. & Dorans, N. J. (1982). The Polyserial Correlation Coefficient. Psychometika, 47, 337-347. }
\author{ John T. Willse }
\examples{
x <- rnorm(500,50,5)
y <- x + rnorm(500,0,2)
x <- x + rnorm(500,0,2)
cor(x,y)

y <- ifelse(y>50,1,0)
cor(x,y)

polyserial(x,y)
polyserial(x,y,ml=FALSE)
}
%\keyword{misc}
