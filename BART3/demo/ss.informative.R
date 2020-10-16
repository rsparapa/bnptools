
library(BART3)

f = function(x)
    5+10*sin(pi*x[ , 1]*x[ , 2]) +  3*x[ , 3]^3
##    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 2000
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
file.='post/ss.cont.bart.rds'
if(!file.exists(file.)) {
    post = mc.gbart(x.train, y.train, x.test,
                    ##sparse=TRUE,
                    mc.cores=B, seed=12)
    saveRDS(post, file.)
} else {
    post = readRDS(file.)
}

shards=8

file.='post/ss.cont.lisa.rds'
if(!file.exists(file.)) {
post2 = ml.gbart(x.train, y.train, x.test, shards=shards,
                 ##sparse=TRUE,
                 mc.cores=B, seed=12)
    saveRDS(post2, file.)
} else {
    post2 = readRDS(file.)
}

post3 = ss.gbart(x.train, y.train, x.test, shards=shards,
                 ##sparse=TRUE,
                 RDSfile='post/ss.cont.gbart',
                 debug=TRUE, mc.cores=B, seed=12)

post$yhat.test.025=apply(post$yhat.test, 2, quantile, probs=0.025)
post$yhat.test.975=apply(post$yhat.test, 2, quantile, probs=0.975)
post2$yhat.test.025=apply(post2$yhat.test, 2, quantile, probs=0.025)
post2$yhat.test.975=apply(post2$yhat.test, 2, quantile, probs=0.975)

for(i in shards) {
    pred=predict(post3[[i]], x.test)
    post3[[i]]$yhat.test.025=apply(pred, 2, quantile, probs=0.025)
    post3[[i]]$yhat.test.975=apply(pred, 2, quantile, probs=0.975)
    post3[[i]]$yhat.test.mean=apply(pred, 2, mean)
}

plot(post$sigma[ , 1], type='l', ylim=c(0, 10),
     ylab=expression(sigma), sub=paste0('N=', N))
for(i in 2:8) lines(post$sigma[ , i])
for(i in 1:shards) lines(post2$sigma[ , i]/sqrt(shards), col=2, lty=2)
##for(i in 1:shards) lines(post2$sigma[ , i], col=2)
##for(i in 1:shards) lines(post3[[shards]]$sigma[ , i]/sqrt(shards), col=8, lty=2)
for(i in 1:shards) lines(post3[[i]]$sigma[ , 1], col=i, lty=3)
for(i in 1:B) lines(post3[[shards]]$sigma[ , i], col=8)
abline(v=100, h=0)
abline(h=1, col='blue')
legend('topright', c('True', 'BART', 'Modified LISA', 'Sequential shards'),
       lty=1, col=c(4, 1:2, 8))
dev.copy2pdf(file='ss-informative1.pdf')

plot(x, f(x.test), type='l', col=4, lwd=2, sub=paste0('N=', N),
     xlab=expression(x[3]), ylab=expression(f(x[3])))
lines(x, post$yhat.test.mean, lwd=2, lty=1)
lines(x, post$yhat.test.025, lwd=2, lty=1)
lines(x, post$yhat.test.975, lwd=2, lty=1)
lines(x, post2$yhat.test.mean, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.025, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.975, lwd=2, lty=2, col=2)
for(i in shards) {
    lines(x, post3[[i]]$yhat.test.mean, lwd=2, lty=3, col=i)
    lines(x, post3[[i]]$yhat.test.025, lwd=2, lty=3, col=i)
    lines(x, post3[[i]]$yhat.test.975, lwd=2, lty=3, col=i)
}
##abline(h=c(0, 5, 10))
##abline(v=0)
legend('topleft', col=c(4, 1:2, 1:8), lty=1,
       legend=c('True', 'BART', 'Modified LISA', paste0(1:8)), lwd=2)
dev.copy2pdf(file='ss-informative2.pdf')

plot(density(c(post3[[shards]]$sigma[-(1:100), ])), ylab='pdf', ylim=c(0, 10),
     xlab=expression(sigma), type='l', xlim=c(0, 3), col=8,
     sub=paste0('N=', N), main='')
##for(i in 1:7) lines(density(c(post3[[i]]$sigma[-(1:100), ])), col=i)
lines(density(c(post2$sigma[-(1:100), ])/sqrt(shards)), col=2, lty=2)
lines(density(c(post$sigma[-(1:100), ])), lty=2)
abline(v=1, h=0, lwd=2, lty=3)
##legend('topleft', col=1:8, lty=1, legend=paste0(1:8), lwd=2)
dev.copy2pdf(file='ss-informative3.pdf')

str(post3$strata)

W=1
w=1
for(i in 2:shards) {
    j=i-1
    W[i]=sum(post3$strata==j)
    w[i]=mean(apply(post3[[i]]$yhat.test, 2, sd)/post3[[i]]$sigma.mean)
}
W=cbind(W, cumsum(W)-1)
W[1, 2]=1
W=cbind(W, sqrt(W[ , 1]/W[ , 2]), sqrt(W[ , 2]/W[ , 1]), w)
W

for(i in 2:shards) {
    j=i-1
    w[i]=mean(apply(post3[[i]]$yhat.test, 2, sd)/
              (post3[[i]]$sigma.mean*W[i, 4]))
}
W=cbind(W, w)
W
