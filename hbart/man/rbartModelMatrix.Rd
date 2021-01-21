\name{rbartModelMatrix}
\alias{rbartModelMatrix}
\docType{package}
\title{Model Matrix for BART}
\description{
rbartModelMatrix takes a data frame of explanatory variables (x)
and turns it into a numeric matrix suitable for BART.
This is used when some of the x variables are factors.
The returned matrix will first have columns for all the numeric variables
in the data frame  and then columns for all the factors expanded into binary
dummy variables.
Note that if a factor has k levels, then BART wants k dummies
(not k-1 as in linear regression).
So for example, a factor with two levels named xf, will result in 
two columns named xf1 and xf2.
}
\usage{
rbartModelMatrix(xdf)
}
\arguments{
\item{xdf}{
Data frame to be turned into a matrix for BART.
}
}
\value{
Matrix from of explanatory variables in the data frame xdf.
}
\examples{
set.seed(99)
xdf = data.frame(x1=1:15,x2=as.factor(c(rep(1,5),rep(2,5),rep(3,5))),
           x3=runif(15),x4=as.factor(c(rep(1,5),rep(2,10))))
print(head(xdf))

xm = rbartModelMatrix(xdf)
print(head(xm))
}
