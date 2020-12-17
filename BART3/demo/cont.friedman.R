
library(BART3)

f = function(x)
    5 + 10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-0.5)^3 + 10*x[ , 4] + 5*x[ , 5]

N = 15000
Q = 5000
sigma = 3 ##y = f(x) + z where z~N(0, sd=sigma)
P = 10    ##number of covariates
B=8
shards=30
ntree=50
nskip=1000

set.seed(21)
x.train=matrix(runif(N*P), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:P)
X.test=matrix(runif(Q*P), Q, P)
dimnames(X.test)[[2]] <- paste0('x', 1:P)

y.train=rnorm(N, f(x.train), sigma)
y.test =rnorm(Q, f(X.test), sigma)

H=100
x=seq(0, 1, length.out=H+1)[-(H+1)]
x.test=matrix(0.5, nrow=H, ncol=P)
x.test[ , 3]=x

file.='post/cont.friedman.rds'
if(!file.exists(file.)) {
    post = mc.gbart(x.train, y.train, x.test,
                    ntree=ntree, nskip=nskip,
                    ##sparse=TRUE,
                    mc.cores=B, seed=12)
    saveRDS(post, file.)
} else {
    post = readRDS(file.)
}

file.='post/cont.lisa.rds'
if(!file.exists(file.)) {
post2 = ml.gbart(x.train, y.train, x.test, shards=shards,
                 ntree=ntree, nskip=nskip,
                 ##sparse=TRUE,
                 mc.cores=B, seed=12)
    saveRDS(post2, file.)
} else {
    post2 = readRDS(file.)
}

file.='post/cont.meta.rds'
if(!file.exists(file.)) {
post5 = ml.gbart(x.train, y.train, x.test, shards=shards,
                 ntree=ntree, nskip=nskip, meta=TRUE,
                 ##sparse=TRUE,
                 mc.cores=B, seed=12)
    saveRDS(post5, file.)
} else {
    post5 = readRDS(file.)
}

file.='post/cont.fweight.rds'
if(!file.exists(file.)) {
post3 = ss.gbart(x.train, y.train, x.test, shards=shards,
                 ntree=ntree, nskip=nskip,
                 ##sparse=TRUE,
                 ##RDSfile='post/ss.cont.gbart',
                 debug=TRUE, mc.cores=B, seed=12)
    saveRDS(post3, file.)
} else {
    post3 = readRDS(file.)
}

file.='post/cont.cweight.rds'
if(!file.exists(file.)) {
post4 = ss.gbart(x.train, y.train, x.test, shards=shards,
                 ntree=ntree, nskip=nskip,
                 ##sparse=TRUE,
                 ##RDSfile='post/ss.cont.gbart',
                 cum.weight=FALSE,
                 debug=TRUE, mc.cores=B, seed=12)
    saveRDS(post4, file.)
} else {
    post4 = readRDS(file.)
}

post$yhat.test.025=apply(post$yhat.test, 2, quantile, probs=0.025)
post$yhat.test.975=apply(post$yhat.test, 2, quantile, probs=0.975)
post2$yhat.test.025=apply(post2$yhat.test, 2, quantile, probs=0.025)
post2$yhat.test.975=apply(post2$yhat.test, 2, quantile, probs=0.975)
post5$yhat.test.025=post5$yhat.test.mean-1.96*post5$yhat.test.sd
post5$yhat.test.975=post5$yhat.test.mean+1.96*post5$yhat.test.sd

for(i in shards) {
    if(i==shards) pred=post3[[i]]$yhat.test
    else pred=predict(post3[[i]], x.test, mc.cores=B)
    post3[[i]]$yhat.test.025=apply(pred, 2, quantile, probs=0.025)
    post3[[i]]$yhat.test.975=apply(pred, 2, quantile, probs=0.975)
    post3[[i]]$yhat.test.mean=apply(pred, 2, mean)
    if(i==shards) pred=post4[[i]]$yhat.test
    else pred=predict(post4[[i]], x.test, mc.cores=B)
    post4[[i]]$yhat.test.025=apply(pred, 2, quantile, probs=0.025)
    post4[[i]]$yhat.test.975=apply(pred, 2, quantile, probs=0.975)
    post4[[i]]$yhat.test.mean=apply(pred, 2, mean)
}

plot(post$sigma[ , 1], type='l', ylim=c(0, 6),
     ylab=expression(sigma), sub=paste0('N=', N))
for(i in 2:8) lines(post$sigma[ , i])
for(i in 1:shards) lines(post2$sigma[ , i]/sqrt(shards), col=2, lty=2)
##for(i in 1:shards) lines(post2$sigma[ , i], col=2)
##for(i in 1:shards) lines(post3[[shards]]$sigma[ , i]/sqrt(shards), col=8, lty=2)
##for(i in 1:shards) lines(post3[[i]]$sigma[ , 1], col=i, lty=3)
for(i in 1:B) lines(post3[[shards]]$sigma[ , i]/sqrt(shards), col=8)
for(i in 1:B) lines(post4[[shards]]$sigma[ , i]/sqrt(shards), col=3)
abline(v=nskip, h=0)
abline(h=sigma, col='blue')
legend('topleft', c('True', 'BART', 'Modified LISA', 'f weight', 'cum. weight'),
       lty=1, col=c(4, 1:3, 8))
dev.copy2pdf(file='cont-friedman-1.pdf')

plot(x, f(x.test), type='l', col=4, lwd=2, sub=paste0('N=', N),
     xlab=expression(x[3]), ylab=expression(f(x[3])))
##     ylim=c(14, 25.5))
##lines(x, f(x.test)-1.96*sigma, type='l', col=4, lwd=2)
##lines(x, f(x.test)+1.96*sigma, type='l', col=4, lwd=2)
lines(x, post$yhat.test.mean, lwd=2, lty=1)
lines(x, post$yhat.test.025, lwd=2, lty=1)
lines(x, post$yhat.test.975, lwd=2, lty=1)
lines(x, post2$yhat.test.mean, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.025, lwd=2, lty=2, col=2)
lines(x, post2$yhat.test.975, lwd=2, lty=2, col=2)
lines(x, post5$yhat.test.mean, lwd=2, lty=2, col=5)
lines(x, post5$yhat.test.025, lwd=2, lty=3, col=5)
lines(x, post5$yhat.test.975, lwd=2, lty=3, col=5)
for(i in shards) {
    lines(x, post3[[i]]$yhat.test.mean, lwd=2, lty=3, col=8)
    lines(x, post3[[i]]$yhat.test.025, lwd=2, lty=3, col=8)
    lines(x, post3[[i]]$yhat.test.975, lwd=2, lty=3, col=8)
    lines(x, post4[[i]]$yhat.test.mean, lwd=2, lty=3, col=3)
    lines(x, post4[[i]]$yhat.test.025, lwd=2, lty=3, col=3)
    lines(x, post4[[i]]$yhat.test.975, lwd=2, lty=3, col=3)
}
##abline(h=c(0, 5, 10))
##abline(v=0)
legend('bottomright', col=c(4, 1:3, 8, 5), lty=1,
       legend=c('True', 'BART', 'Modified LISA', 'f weight', 'cum. weight', 'Meta-analysis'), lwd=2)
dev.copy2pdf(file='cont-friedman-2.pdf')

plot(density(c(post3[[shards]]$sigma[-(1:100), ])/sqrt(shards)), ylab='pdf', ylim=c(0, 10),
     xlab=expression(sigma), type='l', xlim=c(0, 4), col=8,
     sub=paste0('N=', N), main='')
##for(i in 1:7) lines(density(c(post3[[i]]$sigma[-(1:100), ])), col=i)
lines(density(c(post4[[shards]]$sigma[-(1:100), ])/sqrt(shards)), col=3, lty=2)
lines(density(c(post2$sigma[-(1:100), ])/sqrt(shards)), col=2, lty=2)
lines(density(c(post$sigma[-(1:100), ])), lty=2)
abline(v=sigma, h=0, lwd=2, lty=3)
##legend('topleft', col=1:8, lty=1, legend=paste0(1:8), lwd=2)
dev.copy2pdf(file='cont-friedman-3.pdf')

