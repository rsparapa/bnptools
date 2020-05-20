
library(BART3)

f = function(x)
    2.5*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]
N = 10000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates
V = diag(P)
V[3, 4] = 0.8
V[4, 3] = 0.8
L <- chol(V)
set.seed(12)
x.train=matrix(rnorm(N*P), N, P) %*% L
dimnames(x.train)[[2]] <- paste0('x', 1:P)
round(cor(x.train), digits=2)
y.train=(f(x.train)+sigma*rnorm(N)>0)*1
table(y.train)

H=100
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

plot(x.train[ , 3], y.train)
abline(v=0)

B=8
post = mc.gbart(x.train, y.train, x.test, type='pbart',
                sparse=TRUE, mc.cores=B, seed=12)
# set.seed(21)
## post = mc.gbart(x.train, y.train, sparse=TRUE)

post2 = ml.gbart(x.train, y.train, x.test, type='pbart', shards=10,
                 sparse=TRUE, mc.cores=B, seed=12)

post$prob.test.025=apply(post$prob.test, 2, quantile, probs=0.025)
post$prob.test.975=apply(post$prob.test, 2, quantile, probs=0.975)
post2$prob.test.025=apply(post2$prob.test, 2, quantile, probs=0.025)
post2$prob.test.975=apply(post2$prob.test, 2, quantile, probs=0.975)

##pdf('modlisa.pdf')
plot(x, pnorm(f(x.test)),
     type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2, ylim=0:1)
lines(x, post$prob.test.mean, lwd=2, lty=1)
lines(x, post$prob.test.025, lwd=2, lty=2)
lines(x, post$prob.test.975, lwd=2, lty=2)
lines(x, post2$prob.test.mean, lwd=2, col=2)
lines(x, post2$prob.test.025, lwd=2, lty=2, col=2)
lines(x, post2$prob.test.975, lwd=2, lty=2, col=2)
abline(v=0)
legend('topleft', col=c(4, 1, 2), lty=1,
       legend=c('True', '1 shard', '10 shards'), lwd=2)
dev.copy2pdf(file='modlisa-probit.pdf')
##dev.off()
