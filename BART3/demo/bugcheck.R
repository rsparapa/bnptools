
library(BART3)

f = function(x)
    5+10*sin(pi*x[ , 1]*x[ , 2]) +  3*x[ , 3]^3

N = 1000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:P)

y.train=(f(x.train)+sigma*rnorm(N))

H=100
x=seq(-2, 2, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

B=8
post = mc.gbart(x.train, y.train, x.test, mc.cores=B, seed=21)

post2 = mc.wbart(x.train, y.train, x.test, mc.cores=B, seed=21)

set.seed(21)
post3 = gbart(x.train, y.train, x.test)

post$yhat.test.025=apply(post$yhat.test, 2, quantile, probs=0.025)
post$yhat.test.975=apply(post$yhat.test, 2, quantile, probs=0.975)
post2$yhat.test.025=apply(post2$yhat.test, 2, quantile, probs=0.025)
post2$yhat.test.975=apply(post2$yhat.test, 2, quantile, probs=0.975)
post3$yhat.test.025=apply(post3$yhat.test, 2, quantile, probs=0.025)
post3$yhat.test.975=apply(post3$yhat.test, 2, quantile, probs=0.975)

plot(x, f(x.test), type='l', col=4, lwd=2, sub='N=1000',
     xlab=expression(x[3]), ylab=expression(f(x[3])))
lines(x, post$yhat.test.mean, lwd=2, lty=1)
lines(x, post$yhat.test.025, lwd=2, lty=1)
lines(x, post$yhat.test.975, lwd=2, lty=1)
lines(x, post2$yhat.test.mean, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.025, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.975, lwd=2, lty=2, col=2)
lines(x, post3$yhat.test.mean, lwd=2, lty=2, col=3)
lines(x, post3$yhat.test.025, lwd=2, lty=2, col=3)
lines(x, post3$yhat.test.975, lwd=2, lty=2, col=3)
legend('topleft', col=c(4, 1:3), lty=1,
       legend=c('True', 'mc.gbart', 'mc.wbart', 'gbart'), lwd=2)
##dev.copy2pdf(file='ss-informative2.pdf')

post$offset
post3$offset
