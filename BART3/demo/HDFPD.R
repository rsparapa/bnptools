
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

H=20 ### BEWARE H settings and NOT H imputations
x.test=matrix(0, nrow=H, ncol=P)
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test[ , 3]=x

## FPD: no hot-decking
proc.time.=proc.time()
yhat.test=FPD(post, x.train, x.test, 3, mc.cores=B)
print(proc.time()-proc.time.)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)

## HDFPD: naive hot-decking variance
proc.time.=proc.time()
naive=HDFPD(post, x.train, x.test, 3, mult.impute=20, mc.cores=B)
print(proc.time()-proc.time.)

## HDFPD: adjusted hot-decking variance
proc.time.=proc.time()
pred20=HDFPD(post, x.train, x.test, 3, mult.impute=20, hotd.var=TRUE,
             mc.cores=B)
print(proc.time()-proc.time.)

proc.time.=proc.time()
pred5=HDFPD(post, x.train, x.test, 3, mult.impute=5, hotd.var=TRUE, mc.cores=B)
print(proc.time()-proc.time.)

proc.time.=proc.time()
pred2=HDFPD(post, x.train, x.test, 3, mult.impute=2, hotd.var=TRUE, mc.cores=B)
print(proc.time()-proc.time.)

par(mfrow=c(2, 1))
for(i in 1:2) {
    if(i==1) {
        plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', lwd=2,
             xlim=c(-0.5, 0.5), ylim=c(-15, 15))
        } else {
            plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', lwd=2,
                 xlim=c(0.5, 2.5), ylim=c(25, 50))
        legend('topleft', col=0:6, lty=0:6, lwd=2,
           legend=c('Method 1', 'True', 'FPD', 'Naive 20', '20', '5', '2'))
        }
    lines(x, yhat.test.mean, col=2, lwd=2)
    lines(x, yhat.test.025, col=2, lty=2, lwd=2)
    lines(x, yhat.test.975, col=2, lty=2, lwd=2)
    lines(x, naive$yhat.test.lower, col=3, lty=3, lwd=2)
    lines(x, naive$yhat.test.upper, col=3, lty=3, lwd=2)
    lines(x, pred20$yhat.test.mean, col=4, lwd=2)
    lines(x, pred20$yhat.test.lower., col=4, lty=4, lwd=2)
    lines(x, pred20$yhat.test.upper., col=4, lty=4, lwd=2)
    lines(x, pred5$yhat.test.mean, col=5, lwd=2)
    lines(x, pred5$yhat.test.lower., col=5, lty=5, lwd=2)
    lines(x, pred5$yhat.test.upper., col=5, lty=5, lwd=2)
    lines(x, pred2$yhat.test.mean, col=6, lwd=2)
    lines(x, pred2$yhat.test.lower., col=6, lty=6, lwd=2)
    lines(x, pred2$yhat.test.upper., col=6, lty=6, lwd=2)
}
par(mfrow=c(1, 1))
dev.copy2pdf(file='HDFPD-1.pdf')

par(mfrow=c(2, 1))
for(i in 1:2) {
    if(i==1) {
        plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', lwd=2,
             xlim=c(-0.5, 0.5), ylim=c(-15, 15))
        } else {
            plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', lwd=2,
                 xlim=c(0.5, 2.5), ylim=c(25, 50))
        legend('topleft', col=0:6, lty=0:6, lwd=2,
           legend=c('Method 2', 'True', 'FPD', 'Naive 20', '20', '5', '2'))
        }
    lines(x, yhat.test.mean, col=2, lwd=2)
    lines(x, yhat.test.025, col=2, lty=2, lwd=2)
    lines(x, yhat.test.975, col=2, lty=2, lwd=2)
    lines(x, naive$yhat.test.lower, col=3, lty=3, lwd=2)
    lines(x, naive$yhat.test.upper, col=3, lty=3, lwd=2)
    lines(x, pred20$yhat.test.mean, col=4, lwd=2)
    lines(x, pred20$yhat.test.lower., col=4, lty=4, lwd=2)
    lines(x, pred20$yhat.test.upper., col=4, lty=4, lwd=2)
    lines(x, pred5$yhat.test.mean, col=5, lwd=2)
    lines(x, pred5$yhat.test.lower., col=5, lty=5, lwd=2)
    lines(x, pred5$yhat.test.upper., col=5, lty=5, lwd=2)
    lines(x, pred2$yhat.test.mean, col=6, lwd=2)
    lines(x, pred2$yhat.test.lower., col=6, lty=6, lwd=2)
    lines(x, pred2$yhat.test.upper., col=6, lty=6, lwd=2)
}
par(mfrow=c(1, 1))
dev.copy2pdf(file='HDFPD-2.pdf')
