
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

H=20
x.test=matrix(0, nrow=H, ncol=P)
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test[ , 3]=x
yhat.test=FPD(post, x.train, x.test, 3)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)
yhat.test2=HDFPD(post, x.train, x.test, 3, mult.impute=1)
yhat.test2.mean=apply(yhat.test2, 2, mean)
yhat.test2.025=apply(yhat.test2, 2, quantile, probs=0.025)
yhat.test2.975=apply(yhat.test2, 2, quantile, probs=0.975)
yhat.test3=HDFPD(post, x.train, x.test, 3, mult.impute=10)
yhat.test3.mean=apply(yhat.test3, 2, mean)
yhat.test3.025=apply(yhat.test3, 2, quantile, probs=0.025)
yhat.test3.975=apply(yhat.test3, 2, quantile, probs=0.975)

pdf('HDFPD.pdf')
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
lines(x, yhat.test.mean, col=1, lty=1, lwd=2)
lines(x, yhat.test.025, col=1, lty=3, lwd=2)
lines(x, yhat.test.975, col=1, lty=3, lwd=2)
lines(x, yhat.test2.mean/20, col=2, lty=1, lwd=2)
lines(x, yhat.test2.025/20, col=2, lty=3, lwd=2)
lines(x, yhat.test2.975/20, col=2, lty=3, lwd=2)
lines(x, yhat.test3.mean/200, col=3, lty=1, lwd=2)
lines(x, yhat.test3.025/200, col=3, lty=3, lwd=2)
lines(x, yhat.test3.975/200, col=3, lty=3, lwd=2)
legend('topleft', col=c(4, 1, 2, 3), lty=1,
       legend=c('True', 'FPD', '1', '10'), lwd=2)
dev.off()
