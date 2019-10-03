
library(DPM)

N=100
P=10
set.seed(21)
x=c(matrix(runif(N/P), nrow=N/P, ncol=P))
A=quantile(x, (19:21)/100)  ## R
B=.quantile(x, (19:21)/100) ## Rcpp
C=Quantile(x, (19:21)/100)  ## quantile.default
all(A==B)
all(A==C)
all(B==C)
cbind(A, B, C)

all(floor(x)==.floor(x))
all(ceiling(x)==.ceiling(x))
all(sort(x)==.sort(x))
