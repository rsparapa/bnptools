
library(BART3)

f = function(x)
    10*sin(pi*x[ , 1]*x[ , 2]) +  20*x[ , 3]
##    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 1000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates
## P = 10

V = diag(P)
V[3, 4] = 0.8
V[4, 3] = 0.8
## V[5, 6] = 0.8
## V[6, 5] = 0.8
L <- chol(V)
set.seed(12)
x.train=matrix(rnorm(N*P), N, P) %*% L
dimnames(x.train)[[2]] <- paste0('x', 1:P)
round(cor(x.train), digits=2)

y.train=(f(x.train)+sigma*rnorm(N))

B=8
post = mc.gbart(x.train, y.train, sparse=TRUE, mc.cores=B, seed=12)
# set.seed(21)
## post = mc.gbart(x.train, y.train, sparse=TRUE)
sort(post$varprob.mean*P, TRUE)

a=proc.time()
x.test=x.train
yhat.test=SHAP(post, x.train, x.test, 3)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)
print((proc.time()-a)/60)

pdf('SHAP.pdf')
x=x.train[ , 3]
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
points(x, yhat.test.mean, col=1, pch='.')
points(x, yhat.test.025, col=1, pch='.')
points(x, yhat.test.975, col=1, pch='.')
dev.off()
