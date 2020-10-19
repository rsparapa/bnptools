
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

y.train=(1*(f(x.train)+sigma*rnorm(N)>0))

H=100
x=seq(-2, 2, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

B=8
file.='post/ss.probit.bart.rds'
if(!file.exists(file.)) {
    post = mc.gbart(x.train, y.train, x.test, type='pbart',
                    ##sparse=TRUE,
                    mc.cores=B, seed=12)
    saveRDS(post, file.)
} else {
    post = readRDS(file.)
}

shards=8

post3 = ss.zbart(x.train, y.train, x.test, shards=shards, type='pbart',
                 ##sparse=TRUE,
                 ##RDSfile='post/sz.probit.gbart',
                 debug=TRUE, mc.cores=B, seed=12)

table(post3[[shards]]$zdraw)

post$prob.test.025=apply(post$prob.test, 2, quantile, probs=0.025)
post$prob.test.975=apply(post$prob.test, 2, quantile, probs=0.975)

for(i in shards) {
    ##pred=predict(post3[[i]], x.test)
    post3[[i]]$prob.test.025=apply(post3[[i]]$prob.test, 2, quantile, probs=0.025)
    post3[[i]]$prob.test.975=apply(post3[[i]]$prob.test, 2, quantile, probs=0.975)
}

plot(x, pnorm(f(x.test)), type='l', col=4, lwd=2,
     sub=paste0('N=', N),
     ylim=0:1, xlab=expression(x[3]), ylab=expression(f(x[3])))
lines(x, post$prob.test.mean, lwd=2, lty=1)
lines(x, post$prob.test.025, lwd=2, lty=1)
lines(x, post$prob.test.975, lwd=2, lty=1)
for(i in shards) {
    lines(x, post3[[i]]$prob.test.mean, lwd=2, lty=3, col=i)
    lines(x, post3[[i]]$prob.test.025, lwd=2, lty=3, col=i)
    lines(x, post3[[i]]$prob.test.975, lwd=2, lty=3, col=i)
}
##abline(v=0)
legend('topleft', col=c(4, 1, 8), lty=1,
       legend=c('True', 'BART', paste0(8)), lwd=2)
dev.copy2pdf(file='sz-informative-probit.pdf')

