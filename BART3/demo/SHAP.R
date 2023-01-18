
library(BART3)

f = function(x)
    5+10*sin(pi*x[ , 1]*x[ , 2]) +  20*x[ , 3]
##    10*sin(pi*x[ , 1]*x[ , 2]) + 5*x[ , 3]*x[ , 4]^2 + 20*x[ , 5]

N = 500
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
x = x.train[ , 3]
x.train=x.train[order(x), ]
x = x.train[ , 3]
y.train=(f(x.train)+sigma*rnorm(N))

B=8
post = mc.gbart(x.train, y.train, sparse=TRUE, mc.cores=B, seed=12)
# set.seed(21)
## post = mc.gbart(x.train, y.train, sparse=TRUE)
sort(post$varprob.mean*P, TRUE)

## x.test = x.train
H=50
x.test=matrix(0, nrow=H, ncol=P)
x=seq(-3, 3, length.out=H+1)[-(H+1)]
x.test[ , 3]=x

## FPD: no kernel sampling
proc.time.=proc.time()
yhat.test=FPD(post, x.train, x.test, 3, mc.cores=B)
yhat.test=yhat.test-post$offset
print(proc.time()-proc.time.)
yhat.test.mean=apply(yhat.test, 2, mean)
yhat.test.lower=apply(yhat.test, 2, quantile, probs=0.025)
yhat.test.upper=apply(yhat.test, 2, quantile, probs=0.975)

## SHAP: naive kernel sampling variance
file.='SHAP-naive.rds'
if(file.exists(file.)) { naive=readRDS(file.)
} else {
    proc.time.=proc.time()
    naive=SHAP(post, x.train, x.test, 3, mult.impute=50, mc.cores=B)
    print(proc.time()-proc.time.)
    saveRDS(naive, file.)
}

naive=list(SHAP=naive)
naive$yhat.test.mean =apply(naive$SHAP, 2, mean)
naive$yhat.test.lower=apply(naive$SHAP, 2, quantile, probs=0.025)
naive$yhat.test.upper=apply(naive$SHAP, 2, quantile, probs=0.975)
                            
## SHAP: adjusted kernel sampling variance
file.='SHAP-adjust.rds'
if(file.exists(file.)) { adjust=readRDS(file.)
} else {
    proc.time.=proc.time()
    adjust=SHAP(post, x.train, x.test, 3, mult.impute=50, kern.var=TRUE,
                mc.cores=B)
    print(proc.time()-proc.time.)
    saveRDS(adjust, file.)
}

adjust=list(SHAP=adjust)
adjust$yhat.test.mean =apply(adjust$SHAP, 2, mean)
adjust$yhat.test.lower=apply(adjust$SHAP, 2, quantile, probs=0.025)
adjust$yhat.test.upper=apply(adjust$SHAP, 2, quantile, probs=0.975)

pdf(file='SHAP.pdf')
plot(x, 20*x, type='l', xlab='x3', ylab='f(x3)', lwd=2,
     xlim=c(-0.5, 0.5), ylim=c(-15, 15))
lines(x, yhat.test.lower, col=2, lty=2, lwd=2)
lines(x, yhat.test.upper, col=2, lty=2, lwd=2)
lines(x, yhat.test.mean, col=2, lty=2, lwd=2)
lines(x, naive$yhat.test.lower, col=3, lty=3, lwd=2)
lines(x, naive$yhat.test.upper, col=3, lty=3, lwd=2)
lines(x, naive$yhat.test.mean, col=3, lty=3, lwd=2)
lines(x, adjust$yhat.test.lower, col=4, lty=4, lwd=2)
lines(x, adjust$yhat.test.upper, col=4, lty=4, lwd=2)
lines(x, adjust$yhat.test.mean, col=4, lty=4, lwd=2)
legend('topleft', col=0:4, lty=0:4, lwd=2,
       legend=c('Marginal Effects', 'True', 'FPD', 'Naive', 'Adjusted'))
dev.off()

