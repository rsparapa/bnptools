
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
## with forking: Linux/macOS
post = mc.gbart(x.train, y.train, sparse=TRUE, mc.cores=B, seed=12)
## without forking: Windows
## set.seed(21)
## post = gbart(x.train, y.train, sparse=TRUE)

sort(post$varprob.mean*P, TRUE)

H=20
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(x, nrow=H, ncol=P)
yhat.test=FPD(post, x.train=x.train, x.test=x.test, S=3)
yhat.test.mean=apply(yhat.test, 2, mean)

yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)

pdf('FPDex.pdf')
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
lines(x, yhat.test.mean, col=2, lty=2, lwd=2)
lines(x, yhat.test.025, col=2, lty=3, lwd=2)
lines(x, yhat.test.975, col=2, lty=3, lwd=2)
legend('topleft', col=c(4, 2, 1), lty=c(1, 2),
       legend=c('True', 'FPD'), lwd=2)
dev.off()

marg=FPDK(post, x.train=x.train, x.test=x.test, S=3, mult.impute=20)
yhat.test.meanK=apply(marg$yhat.test, 2, mean)
yhat.test.025K=apply(marg$yhat.test, 2, quantile, probs=0.025)
yhat.test.975K=apply(marg$yhat.test, 2, quantile, probs=0.975)

pdf('FPDKex.pdf')
plot(x, 20*x, type='n', xlab='x3', ylab='f(x3)', col=4, lwd=2)
lines(x, yhat.test.mean, col=2, lty=2, lwd=2)
lines(x, yhat.test.025, col=2, lty=3, lwd=2)
lines(x, yhat.test.975, col=2, lty=3, lwd=2)
lines(x, yhat.test.meanK, col=4, lty=2, lwd=2)
lines(x, yhat.test.025K, col=4, lty=3, lwd=2)
lines(x, yhat.test.975K, col=4, lty=3, lwd=2)
legend('topleft', col=c(4, 2, 1), lty=c(1, 2),
       legend=c('FPDK', 'FPD'), lwd=2)
dev.off()
