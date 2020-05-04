
library(BART3)

f = function(x) #only the first 5 matter
    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 1000
sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
P = 10       #number of covariates

V = diag(P)
V[5, 6] = 0.8
V[6, 5] = 0.8
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


pdf('shapley.pdf')
H=20
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test=matrix(x, nrow=H, ncol=P)
yhat.test=SHAP(post$treedraws, x.train, x.test, 5)
yhat.test.mean=apply(yhat.test, 2, mean)

## Friedman's partial dependence function
pred=list()
for(h in 1:H) {
    x.test=x.train
    x.test[ , 5]=x[h]
    pred[[h]]=predict(post, x.test, mc.cores=B)
    if(h==1) yhat.test2=cbind(apply(pred[[h]], 1, mean))
    else yhat.test2=cbind(yhat.test2, cbind(apply(pred[[h]], 1, mean)))
}
yhat.test2.mean=apply(yhat.test2, 2, mean)

plot(x, yhat.test.mean, type='l', xlab='x5', ylab='f(x5)', lty=3)
lines(x, yhat.test2.mean, col=4, lty=2)
lines(x, 20*x, col=2)

dev.off()
##dev.copy2pdf(file='shapley.pdf')
