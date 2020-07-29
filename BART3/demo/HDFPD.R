
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
## FPD: truth
yhat.test=FPD(post, x.train, x.test, 3)
yhat.test.mean=apply(yhat.test, 2, mean)
fvar.test=apply(yhat.test, 2, var)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)

## HDFPD: approximation
yhat.test20=HDFPD(post, x.train, x.test, 3, mult.impute=20)
yhat.test20.mean=apply(yhat.test20, 2, mean)
Yhat.test20.mean=matrix(yhat.test20.mean, byrow=TRUE,
                        nrow=nrow(yhat.test20),
                        ncol=ncol(yhat.test20))
yhat.test20.025=apply(yhat.test20, 2, quantile, probs=0.025)
yhat.test20.975=apply(yhat.test20, 2, quantile, probs=1-0.025)

## HDFPD: approximation with correction
pred=HDFPD(post, x.train, x.test, 3, mult.impute=20, hdvar=TRUE)
fvar.test20=apply(pred$yhat.test, 2, var)
hdvar.test20=apply(pred$hdvar.test, 2, mean)
hdvar.test20/H
##c=(1-(hdvar.test20/fvar.test20)/H)^(-0.5)
fvar.test20.=fvar.test20-(hdvar.test20/H)
fsd.test20.=sqrt(fvar.test20.)

yhat.test20.=Yhat.test20.mean+
    (yhat.test20-Yhat.test20.mean)*
    sqrt(fvar.test20./fvar.test20)
yhat.test20.025.=apply(yhat.test20., 2, quantile, probs=0.025)
yhat.test20.975.=apply(yhat.test20., 2, quantile, probs=1-0.025)

##pdf('HDFPD.pdf')
par(mfrow=c(2, 1))
for(i in 1:2) {
    if(i==1) {
        plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2,
             xlim=c(-1, 1), ylim=c(-15, 15))
        legend('topleft', col=c(4, 1, 2, 3), lty=1,
           legend=c('True', 'FPD', 'Method 1', 'no correction'), lwd=2)
        } else
            plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2,
                 xlim=c(1, 3), ylim=c(15, 45))
    lines(x, yhat.test.mean, col=1, lty=2, lwd=2)
    lines(x, yhat.test.025, col=1, lty=3, lwd=2)
    lines(x, yhat.test.975, col=1, lty=3, lwd=2)
    lines(x, yhat.test20.mean, col=2, lty=2, lwd=2)
    lines(x, yhat.test20.025, col=3, lty=3, lwd=2)
    lines(x, yhat.test20.975, col=3, lty=3, lwd=2)
    lines(x, yhat.test20.mean-1.96*fsd.test20., col=2, lty=3, lwd=2)
    lines(x, yhat.test20.mean+1.96*fsd.test20., col=2, lty=3, lwd=2)
}
par(mfrow=c(1, 1))
dev.copy2pdf(file='HDFPD-1.pdf')

par(mfrow=c(2, 1))
for(i in 1:2) {
    if(i==1) {
        plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2,
             xlim=c(-1, 1), ylim=c(-15, 15))
        legend('topleft', col=c(4, 1, 2, 3), lty=1,
           legend=c('True', 'FPD', 'Method 2', 'no correction'), lwd=2)
        } else
            plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2,
                 xlim=c(1, 3), ylim=c(15, 45))
    lines(x, yhat.test.mean, col=1, lty=2, lwd=2)
    lines(x, yhat.test.025, col=1, lty=3, lwd=2)
    lines(x, yhat.test.975, col=1, lty=3, lwd=2)
    lines(x, yhat.test20.mean, col=2, lty=2, lwd=2)
    lines(x, yhat.test20.025, col=3, lty=3, lwd=2)
    lines(x, yhat.test20.975, col=3, lty=3, lwd=2)
    lines(x, yhat.test20.025., col=2, lty=3, lwd=2)
    lines(x, yhat.test20.975., col=2, lty=3, lwd=2)
}
par(mfrow=c(1, 1))
dev.copy2pdf(file='HDFPD-2.pdf')
