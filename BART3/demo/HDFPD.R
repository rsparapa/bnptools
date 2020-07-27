
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
fvar.test=apply(yhat.test, 2, var)
yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)
## yhat.test10=HDFPD(post, x.train, x.test, 3, mult.impute=10)
## yhat.test10.mean=apply(yhat.test10, 2, mean)
## yhat.test10.025=apply(yhat.test10, 2, quantile, probs=0.025)
## yhat.test10.975=apply(yhat.test10, 2, quantile, probs=1-0.025)
## yhat.test20=HDFPD(post, x.train, x.test, 3, mult.impute=20)
## yhat.test20.mean=apply(yhat.test20, 2, mean)
yhat.test20.025=apply(yhat.test20, 2, quantile, probs=0.025)
yhat.test20.975=apply(yhat.test20, 2, quantile, probs=1-0.025)

pred=HDFPD(post, x.train, x.test, 3, mult.impute=20, hdvar=TRUE)
c=(1-(hdvar.test20/fvar.test20)/20)^(-0.5)
fvar.test20=apply(pred$yhat.test, 2, var)
hdvar.test20=apply(pred$hdvar.test, 2, mean)
hdvar.test20

yhat.test20.025.=0
yhat.test20.975.=0
for(i in 1:H) {
yhat.test20.025.[i]=quantile(pred$yhat.test[ , i], probs=0.025*c[i])
yhat.test20.975.[i]=quantile(pred$yhat.test[ , i], probs=1-0.025*c[i])
}

## yhat.test2.025=apply(yhat.test2, 2, quantile, probs=0.025*sqrt(1000/20))
## yhat.test2.975=apply(yhat.test2, 2, quantile, probs=1-0.025*sqrt(1000/20))
## yhat.test3=HDFPD(post, x.train, x.test, 3, mult.impute=10)
## yhat.test3.mean=apply(yhat.test3, 2, mean)
## yhat.test3.025=apply(yhat.test3, 2, quantile, probs=0.025)
## yhat.test3.975=apply(yhat.test3, 2, quantile, probs=0.975)
## yhat.test4=HDFPD(post, x.train, x.test, 3, mult.impute=20)
## yhat.test4.mean=apply(yhat.test4, 2, mean)
## yhat.test4.025=apply(yhat.test4, 2, quantile, probs=0.025)
## yhat.test4.975=apply(yhat.test4, 2, quantile, probs=0.975)
## yhat.test5=HDFPD(post, x.train, x.test, 3, mult.impute=40)
## yhat.test5.mean=apply(yhat.test5, 2, mean)
## yhat.test5.025=apply(yhat.test5, 2, quantile, probs=0.025)
## yhat.test5.975=apply(yhat.test5, 2, quantile, probs=0.975)

##pdf('HDFPD.pdf')
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
lines(x, yhat.test.mean, col=1, lty=2, lwd=2)
lines(x, yhat.test.025, col=1, lty=3, lwd=2)
lines(x, yhat.test.975, col=1, lty=3, lwd=2)
lines(x, yhat.test20.mean, col=2, lty=2, lwd=2)
lines(x, yhat.test20.025., col=2, lty=3, lwd=2)
lines(x, yhat.test20.975., col=2, lty=3, lwd=2)
## lines(x, yhat.test10.mean, col=3, lty=2, lwd=2)
lines(x, yhat.test20.025, col=3, lty=3, lwd=2)
lines(x, yhat.test20.975, col=3, lty=3, lwd=2)
legend('topleft', col=c(4, 1, 2, 3), lty=1,
       legend=c('True', 'FPD', 'correction', 'no correction'), lwd=2)
## lines(x, yhat.test3.mean, col=3, lty=1, lwd=2)
## lines(x, yhat.test3.025, col=3, lty=3, lwd=2)
## lines(x, yhat.test3.975, col=3, lty=3, lwd=2)
## lines(x, yhat.test4.mean, col=5, lty=1, lwd=2)
## lines(x, yhat.test4.025, col=5, lty=3, lwd=2)
## lines(x, yhat.test4.975, col=5, lty=3, lwd=2)
## lines(x, yhat.test5.mean, col=6, lty=1, lwd=2)
## lines(x, yhat.test5.025, col=6, lty=3, lwd=2)
## lines(x, yhat.test5.975, col=6, lty=3, lwd=2)
## legend('topleft', col=c(4, 1, 2, 3, 5), lty=1,
##        legend=c('True', 'FPD', '1', '10', '20', '40'), lwd=2)
##dev.off()
dev.copy2pdf(file='HDFPD.pdf')
