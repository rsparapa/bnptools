
f = function(x)
    -10+10*x[ , 1]+5*x[ , 2]+10*sin(pi*x[ , 4]*x[ , 5])+20*(x[ , 6]-0.5)^3

P = 50
N = 200
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

check = gbart(X., y, impute.mult=1:3)
table(check, miss)
