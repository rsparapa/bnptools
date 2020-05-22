
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
##x.train=x.train[order(x.train[ , 3]), ]
y.train=(f(x.train)+sigma*rnorm(N))

B=8
post = mc.gbart(x.train, y.train, sparse=TRUE, mc.cores=B, seed=12)
# set.seed(21)
## post = mc.gbart(x.train, y.train, sparse=TRUE)
sort(post$varprob.mean*P, TRUE)

pdf('HD.pdf')
for(mc.cores in c(8)) {
    a=proc.time()
    x=x.train[ , 3]
    x.test = x.train
    yhat.test=HD(post, x.train, x.test, 3, mc.cores=mc.cores)
    yhat.test.mean=apply(yhat.test, 2, mean)
    yhat.test.025=apply(yhat.test, 2, quantile, probs=0.025)
    yhat.test.975=apply(yhat.test, 2, quantile, probs=0.975)
    plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
    points(x, yhat.test.mean, col=2, lwd=2, pch='.')
    points(x, yhat.test.025, col=1, lwd=2, pch='.')
    points(x, yhat.test.975, col=1, lwd=2, pch='.')
    legend('topleft', col=c(4, 2, 1), lty=c(1, 2),
           legend=c('True', 'HD'), lwd=2)
    print((proc.time()-a)/60)
}
dev.off()
