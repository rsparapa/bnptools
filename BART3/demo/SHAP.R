
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
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(x, nrow=H, ncol=P)
yhat.test=SHAP(post, x.train, x.test, 3)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)
## yhat.test2=FPD(post, x.train, x.test, 3)
## yhat.test2.mean=apply(yhat.test2, 2, mean)
## yhat.test2.025=apply(yhat.test2, 2, quantile, probs=0.025)
## yhat.test2.975=apply(yhat.test2, 2, quantile, probs=0.975)

yhat.test3=HD(post, x.train, x.test, 3)
yhat.test3.mean=apply(yhat.test3, 2, mean)
yhat.test3.025=apply(yhat.test3, 2, quantile, probs=0.025)
yhat.test3.975=apply(yhat.test3, 2, quantile, probs=0.975)

pdf('SHAP.pdf')
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
lines(x, yhat.test.mean, lwd=2, lty=3)
lines(x, yhat.test.025, lty=2, lwd=2)
lines(x, yhat.test.975, lty=2, lwd=2)
## lines(x, yhat.test2.mean, col=2, lwd=2)
## lines(x, yhat.test2.025, col=2, lty=2, lwd=2)
## lines(x, yhat.test2.975, col=2, lty=2, lwd=2)
lines(x, yhat.test3.mean, col=2, lwd=2, lty=4)
lines(x, yhat.test3.025, col=2, lty=2, lwd=2)
lines(x, yhat.test3.975, col=2, lty=2, lwd=2)
## legend('topleft', col=c(4, 1, 2, 3), lty=c(1, 3, 1, 4),
##        legend=c('True', 'SHAP', 'FPD', 'HD'), lwd=2)
legend('topleft', col=c(4, 1, 2), lty=c(1, 2, 4),
       legend=c('True', 'SHAP', 'HD'), lwd=2)
dev.off()
