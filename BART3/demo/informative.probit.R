
library(BART3)

f = function(x)
    10*sin(pi*x[ , 1]*x[ , 2]) +  3*x[ , 3]^3
##    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 10000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates
## P = 10

V = diag(P)
## V[3, 4] = 0.8
## V[4, 3] = 0.8
L <- chol(V)
set.seed(12)
x.train=matrix(rnorm(N*P), N, P) %*% L
dimnames(x.train)[[2]] <- paste0('x', 1:P)
round(cor(x.train), digits=2)

y.train=((f(x.train)+sigma*rnorm(N))>0)

H=100
x=seq(-2, 2, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

B=8
file.='post/probit.bart.rds'
if(!file.exists(file.)) {
    post = mc.gbart(x.train, y.train, x.test, type='pbart',
                    sparse=TRUE, mc.cores=B, seed=12)
    saveRDS(post, file.)
} else {
    post = readRDS(file.)
}

## file.='post/probit.lisa.rds'
## if(!file.exists(file.)) {
## post2 = ml.gbart(x.train, y.train, x.test, type='pbart', shards=10,
##                  sparse=TRUE, mc.cores=B, seed=12)
##     saveRDS(post2, file.)
## } else {
##     post2 = readRDS(file.)
## }

post3 = ip.gbart(x.train, y.train, x.test, type='pbart', shards=10,
                 sparse=TRUE, mc.cores=B, seed=12)

post$prob.test.025=apply(post$prob.test, 2, quantile, probs=0.025)
post$prob.test.975=apply(post$prob.test, 2, quantile, probs=0.975)
## post2$prob.test.025=apply(post2$prob.test, 2, quantile, probs=0.025)
## post2$prob.test.975=apply(post2$prob.test, 2, quantile, probs=0.975)
post3$prob.test.025=apply(post3$prob.test, 2, quantile, probs=0.025)
post3$prob.test.975=apply(post3$prob.test, 2, quantile, probs=0.975)

## plot(post$sigma[ , 1], type='l', ylim=c(0, 10),
##      ylab=expression(sigma), sub='N=10000')
## for(i in 2:8) lines(post$sigma[ , i])
## for(i in 1:10) lines(post2$sigma[ , i]/sqrt(10), col=2)
## for(i in 1:8) lines(post3$sigma[ , i], col=8)
## abline(v=100, h=0)
## abline(h=1, col='blue')
## legend('topright', c('BART', 'Modified LISA',
##                      'Sequential shards', 'True'),
##        lty=1, col=c(1:2, 8, 4))
plot(x, pnorm(f(x.test)), type='l', col=4, lwd=2, sub='N=10000', ylim=0:1, 
     xlab=expression(x[3]),
     ylab=expression(Phi(f(x[3]))))
lines(x, post$prob.test.mean, lwd=2, lty=1)
lines(x, post$prob.test.025, lwd=2, lty=2)
lines(x, post$prob.test.975, lwd=2, lty=2)
## lines(x, post2$prob.test.mean, lwd=2, col=2)
## lines(x, post2$prob.test.025, lwd=2, lty=2, col=2)
## lines(x, post2$prob.test.975, lwd=2, lty=2, col=2)
lines(x, post3$prob.test.mean, lwd=2, col=8)
lines(x, post3$prob.test.025, lwd=2, lty=2, col=8)
lines(x, post3$prob.test.975, lwd=2, lty=2, col=8)
legend('topleft', col=c(4, 1, 8), lty=1,
       legend=c('True', 'BART', 
                'Sequential shards'), lwd=2)
dev.copy2pdf(file='probit-cube.pdf')

