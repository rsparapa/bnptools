
library(BART3)

N = 1000
P = 5       ##number of covariates
C = 8

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
x.train=x.train[order(x.train[ , 1]), ]
y.true = x.train[ , 1]^3
z = rbinom(N, 1, 0.5)
epsilon = ((-1)^z)*abs(rnorm(N)*z+rnorm(N, 0, 2)*(1-z))
y.train=y.true+epsilon

plot(x.train[ , 1], y.true, type='l')
points(x.train[ , 1], y.train)

##run BART with C cores in parallel
post = mc.gbart(x.train, y.train, mc.cores=C, seed=99)

lines(x.train[ , 1], post$yhat.train.mean, col=2)

set.seed(21)
post2 = liobart(x.train, y.train, nskip=0, ndpost=1)
