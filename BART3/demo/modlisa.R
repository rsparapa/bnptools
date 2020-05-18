
library(BART3)

f = function(x)
    10*sin(pi*x[ , 1]*x[ , 2]) +  20*x[ , 3]
##    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 10000
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

H=20
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

post2 = ml.gbart(x.train, y.train, x.test, shards=10,
                 sparse=TRUE, mc.cores=B, seed=12)

pdf('modlisa.pdf')
plot(post$sigma[ , 1], type='l', ylim=c(0, 10), ylab='SD')
for(i in 1:8) {
    if(i>1) lines(post$sigma[ , i])
    lines(apply(post2$sigma[ , seq(i, 80, 10)], 1, mean)/sqrt(10), col=2)
}
abline(v=100, h=0)
abline(h=1, col='blue')
legend('topright', c('N=10000', '1 shard', '10 shards', 'True'),
       lty=1, col=c(0:2, 4))
dev.off()

## yhat.test=HD(post, x.train, x.test, 3)
## yhat.test.mean=apply(yhat.test, 2, mean)
## yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
## yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)

## pdf('modlisa.pdf')
## plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
## lines(x, yhat.test.mean, lwd=2, lty=3)
## lines(x, yhat.test.025, lty=2, lwd=2)
## lines(x, yhat.test.975, lty=2, lwd=2)
## legend('topleft', col=c(4, 1), lty=c(1, 2),
##        legend=c('True', 'HD'), lwd=2)
## dev.off()
