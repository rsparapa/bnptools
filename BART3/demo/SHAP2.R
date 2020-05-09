
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

H=10
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H^2, ncol=P)
h=0
for(i in 1:H)
    for(j in 1:H) {
        h=h+1
        x.test[h, 1]=x[i]
        x.test[h, 2]=x[j]
}
h
dimnames(x.test)[[2]]=paste0('x', 1:4)        
yhat.test=SHAP2(post, x.train, x.test, 1:2)
yhat.test.mean=apply(yhat.test, 2, mean)

z.true=matrix(nrow=H, ncol=H)
zhat.test=matrix(nrow=H, ncol=H)
h=0
for(i in 1:H)
    for(j in 1:H) {
        h=h+1
        z.true[i, j]=c(f(cbind(x[i], x[j], 0, 0)))
        zhat.test[i, j]=yhat.test.mean[h]
}
h

pdf('SHAP2.pdf')
##par(mfrow=c(2, 1))
contour(x, x, z.true, xlab='x1', ylab='x2', col='blue', nlevels=7)
contour(x, x, zhat.test, add=TRUE, col='red', nlevels=7)
dev.off()
par(mfrow=c(1, 1))

## Friedman's partial dependence function
## pred=list()
## for(h in 1:H) {
##     x.test=x.train
##     x.test[ , 3]=x[h]
##     pred[[h]]=predict(post, x.test, mc.cores=B)
##     if(h==1) yhat.test2=cbind(apply(pred[[h]], 1, mean))
##     else yhat.test2=cbind(yhat.test2, cbind(apply(pred[[h]], 1, mean)))
## }
## yhat.test2.mean=apply(yhat.test2, 2, mean)

## pdf('shapley.pdf')
## plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', col=4, lwd=2)
## lines(x, yhat.test.mean, col=2, lty=3, lwd=2)
## lines(x, yhat.test2.mean, lty=2, lwd=2)
## legend('topleft', col=c(4, 2, 1), lty=c(1, 3, 2),
##        legend=c('True', 'SHAP', 'Friedman'), lwd=2)
## dev.off()
##dev.copy2pdf(file='shapley.pdf')
