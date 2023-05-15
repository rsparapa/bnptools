
library(nftbart)

B=8
options(mc.cores=B)

data(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1

## this study reports time in days rather than weeks or months
times=lung$time
times=times/7  ## weeks

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])
## lung$sex:        Male=1 Female=2

set.seed(99)
post=nft(x.train, times, delta, K=0, mask=0.95)

print(Cindex(-post$f.train.mean, times, delta))

x.test = rbind(x.train, x.train)
x.test[ , 2]=rep(1:2, each=N)
K=100
events=quantile(times[delta==1], (0:(K-1))/(K-1))

a=proc.time()
pred = predict(post, x.test, K=K, events=events, XPtr=TRUE)
print((proc.time()-a)/60)

NK=N*K
RR = matrix(nrow=length(post$s.train.mask), ncol=K)
for(j in 1:K) {
    h=seq(j, NK, K)
    RR[ , j]=apply(pred$haz.test[ , h]/pred$haz.test[ , NK+h], 1, mean)
}

RR.mean = apply(RR, 2, mean)
RR.lower = apply(RR, 2, quantile, 0.025)
RR.upper = apply(RR, 2, quantile, 0.975)
RR.median. = apply(RR, 1, median)
RR.lower. = quantile(RR.median., probs=0.025)
RR.upper. = quantile(RR.median., probs=0.975)
RR.median. = median(RR.median.)

plot(events, RR.mean, type='l', log='y',
     ylim=c(0.02, 50), ylab='RR(t, x)', xlab='t (weeks)',
     sub="Friedman's partial dependence function")
lines(events, RR.lower, lty=2)
lines(events, RR.upper, lty=2)
abline(h=1, col=8)
abline(h=RR.median., col=2)
abline(h=c(RR.lower., RR.upper.), col=2, lty=2)
legend('bottomleft', c('M vs. F', 'Time-varying', 'Proportionality'),
       col=c(0, 1, 2), lty=1)



