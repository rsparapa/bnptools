
library(BART3)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) ## only the first 5 matter
    10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-0.5)^2 + 10*x[ , 4] + 5*x[ , 5]

P = 125       ## number of covariates: signal+noise
C = 60        ## number of continuous: noise
C. = C+5      ## number of continuous: signal+noise
D = P-C.      ## number of dichotomous noise
nskip = 1000  ## burn-in discarded
N = 100       ## sample size
H = 20        ## number of trees
T = 40        ## number of steps

set.seed(34)
x.train=matrix(c(runif(N*C.), rbinom(N*D, 1, 0.5)), N, P)
y.train=rnorm(N, f(x.train))

sparse <- tsvs(x.train, y.train, nskip = nskip, ntree = H, T = T,
               rds.file = 'sparse-tsvs.rds', pdf.file = 'sparse-tsvs.pdf')
