
library(BART3)

f = function(x)
    -15+20*x[ , 1]+10*x[ , 2]+10*sin(pi*x[ , 4]*x[ , 5])+20*(x[ , 6]-0.5)^3

P = 50
N = 100
set.seed(12)
X = matrix(runif(P*N), nrow=N, ncol=P)
x = t(rmultinom(N, 1, c(0.17, 0.49, 0.34)))
X[ , 1:3] = x
y = rnorm(N, f(X))

summary(y[x[ , 1]==1])
summary(y[x[ , 2]==1])
summary(y[x[ , 3]==1])

miss = rbinom(N, 1, 0.4)
x. = x
x.[which(miss==1), 1:3]=NA
X. = X
X.[ , 1:3] = x.

set.seed(34)
post = gbart(X., y, impute.mult=1:3, sparse=TRUE)

post2 = mc.gbart(X., y, sparse=TRUE, seed=22)

post$varprob.mean[1:6]  ## missing imputation
post2$varprob.mean[1:6] ## vs. hot-decking

cor(post$yhat.train.mean, y)^2
cor(post2$yhat.train.mean, y)^2

## set.seed(8)
## post3 = gbart(X., y, sparse=TRUE)
## table(post3[1, ], post3[2, ])
## table(post3[1, ], post3[3, ])
## table(post3[2, ], post3[3, ])
