
library(nftbart)

B=8
options(mc.cores=B)
##B=getOption('mc.cores', 1)
##figures = getOption('figures', default='NONE')

data(lung)
str(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1
table(delta)

## this study reports time in days
times=lung$time
times=times/7  ## weeks
summary(times)

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])
## lung$sex:        Male=1 Female=2

a=proc.time()
file.='lung.rds'
if(file.exists(file.)) {
    post=readRDS(file.)
    XPtr=FALSE
} else {
    set.seed(99)
    post=nft(x.train, times, delta, tc=B, K=0)
    XPtr=TRUE
    ##saveRDS(post, file.)
}
print((proc.time()-a)/60)

pred=predict(post, x.train, tc=B, XPtr=XPtr,
                 seed=21, mask=TRUE)

print(Cindex(pred$logt.test.mean, times, delta))

