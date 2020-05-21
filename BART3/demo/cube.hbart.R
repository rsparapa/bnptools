
library(BART3)

N = 1000
P = 5       #number of covariates
ndpost = 1000
nskip = 100
C = 8

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
x.train=x.train[order(x.train[ , 1]), ]
yhat.train = x.train[ , 1]^3
shat.train = exp(x.train[ , 1])
y.train=rnorm(N, yhat.train, shat.train)

plot(x.train[ , 1], y.train)

post = hbart(x.train, y.train,
                seed=99, ndpost=100, nskip=0)

