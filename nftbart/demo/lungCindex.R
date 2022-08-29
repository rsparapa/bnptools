
library(nftbart)

B=8
options(mc.cores=B)

data(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1

## this study reports time in days
times=lung$time
times=times/7  ## weeks

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])
## lung$sex:        Male=1 Female=2

set.seed(99)
post=nft(x.train, times, delta, K=0)
pred=predict(post, x.train, XPtr=TRUE, seed=21)
print(Cindex(pred$logt.test.mean, times, delta))

