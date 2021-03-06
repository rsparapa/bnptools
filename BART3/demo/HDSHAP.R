
library(BART3)

f=function(x)
    10*(sin(2*pi*x[ , 1]*x[ , 2])+2*x[ , 3]-1)

N = 500
sigma = 1 ##y = f(x) + sigma*z where z~N(0, 1)
##P = 4   ##number of covariates
P = 4
S = 3
C = 2 ## sample C combinations
set.seed(12)
## V = diag(P)
## V[3, 4] = 0.8
## V[4, 3] = 0.8
## ## V[5, 6] = 0.8
## ## V[6, 5] = 0.8
## L <- chol(V)
## x.train=matrix(rnorm(N*P), N, P) %*% L
## dimnames(x.train)[[2]] <- paste0('x', 1:P)
## round(cor(x.train), digits=2)
x.train=matrix(runif(N*P), nrow=N, ncol=P)
x = x.train[ , S]
x.train=x.train[order(x), ]
x = x.train[ , S]
## y.train=(f(x.train)+sigma*rnorm(N))
y.train=f(x.train)

B=8
post = mc.gbart(x.train, y.train, sparse=TRUE, mc.cores=B, seed=12)
# set.seed(21)
## post = mc.gbart(x.train, y.train, sparse=TRUE)
sort(post$varprob.mean*P, TRUE)

x.test = x.train
## H=20
## x.test=matrix(0, nrow=H, ncol=P)
## x=seq(-3, 3, length.out=H+1)[-(H+1)]
## x.test[ , 3]=x

## FPD: no hot-decking
proc.time.=proc.time()
yhat.test=FPD(post, x.train, x.test, S, mc.cores=B)-post$offset
print(proc.time()-proc.time.)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)

## HDSHAP: naive hot-decking variance
proc.time.=proc.time()
naive=HDSHAP(post, x.train, x.test, S, mult.impute=20, mc.cores=B)
print(proc.time()-proc.time.)

## HDSHAP: adjusted hot-decking variance
proc.time.=proc.time()
pred20=HDSHAP(post, x.train, x.test, S, comb.draw=C, mult.impute=20, hotd.var=TRUE, mc.cores=B)
print(proc.time()-proc.time.)

proc.time.=proc.time()
pred5=HDSHAP(post, x.train, x.test, S, comb.draw=C, mult.impute=5, hotd.var=TRUE, mc.cores=B)
print(proc.time()-proc.time.)

proc.time.=proc.time()
pred2=HDSHAP(post, x.train, x.test, S, comb.draw=C, mult.impute=2, hotd.var=TRUE, mc.cores=B)
print(proc.time()-proc.time.)

        plot(x, f(cbind(0, 0, x)), type='l', xlab='x3', ylab='f(x3)', lwd=2,
                     ylim=c(-5, 5), xlim=c(0.2, 0.8))
        legend('topleft', lwd=2,
               col=1:6, lty=1:6,
           legend=c('True', 'FPD', 'Naive 20', '20', '5', '2'))
    lines(x, yhat.test.mean, col=2, lwd=2)
    lines(x, yhat.test.025, col=2, lty=2, lwd=2)
    lines(x, yhat.test.975, col=2, lty=2, lwd=2)
    lines(x, naive$yhat.test.lower, col=3, lty=3, lwd=2)
    lines(x, naive$yhat.test.upper, col=3, lty=3, lwd=2)
    lines(x, pred20$yhat.test.mean, col=4, lwd=2)
    lines(x, pred20$yhat.test.lower, col=4, lty=4, lwd=2)
    lines(x, pred20$yhat.test.upper, col=4, lty=4, lwd=2)
    lines(x, pred5$yhat.test.mean, col=5, lwd=2)
    lines(x, pred5$yhat.test.lower, col=5, lty=5, lwd=2)
    lines(x, pred5$yhat.test.upper, col=5, lty=5, lwd=2)
    lines(x, pred2$yhat.test.mean, col=6, lwd=2)
    lines(x, pred2$yhat.test.lower, col=6, lty=6, lwd=2)
    lines(x, pred2$yhat.test.upper, col=6, lty=6, lwd=2)
dev.copy2pdf(file='HDSHAP.pdf')

