
library(nftbart)

B=8
##B=getOption('mc.cores', 1)
##figures = getOption('figures', default='NONE')

data(lung)
str(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1
table(delta)

## this study reports time in days rather than weeks or months
times=lung$time
times=times/7  ## weeks
summary(times)

## matrix of covariates
x.train=cbind(lung[ , -(1:3)])

D=c(8, 2)
M=1000
T=50
B=8
C=0.5
namesX=dimnames(x.train)[[2]] 
print(namesX)
P=ncol(x.train)
file.='TSVSlung.rds'
if(file.exists(file.)) {
    prob=readRDS(file.)
} else {
    a=matrix(0, nrow=T, ncol=P)
    b=matrix(0, nrow=T, ncol=P)
    S=matrix(0, nrow=T, ncol=P)
    dimnames(S)[[2]]=namesX
    theta=matrix(nrow=T, ncol=P)
    gamma=matrix(0, nrow=T, ncol=P)
    prob=matrix(nrow=T, ncol=P)
    dimnames(prob)[[2]]=namesX
    varcount=matrix(0, nrow=T, ncol=P)
    dimnames(varcount)[[2]]=namesX
    for(t in 1:T) {
        set.seed(T+t)
        print(t)
        for(i in 1:P) {
            if(t==1) {
                a[t, i]=1
                b[t, i]=1
            } else {
                a[t, i]=a[t-1, i]
                b[t, i]=b[t-1, i]
            }
            theta[t, i]=rbeta(1, a[t, i], b[t, i])
            if(theta[t, i]>=C) {
                S[t, i]=1
            }
            ##else S[t, i]=0
        }
print(S[t, ])
##print(S.[t, ])
        set.seed(t)
        x.train.=cbind(x.train[ , S[t, ]==1])
        dimnames(x.train.)[[2]]=namesX[S[t, ]==1]
        post=nft(x.train., times, delta, tc=B, ndpost=M, ntree=D)
        namesV=dimnames(post$f.varcount)[[2]]
        for(i in 1:P) {
            if(S[t, i]==1) {
                k.=which(namesX[i]==namesV)
                j=post$f.varcount[M, k.]+post$s.varcount[M, k.]
                if(j>0) {
                    varcount[t, i]=j
                    gamma[t, i]=1
                    a[t, i]=a[t, i]+1
                } else {
                    b[t, i]=b[t, i]+1
                }
            } else {
                b[t, i]=b[t, i]+1
            }
            prob[t, i]=a[t, i]/(a[t, i]+b[t, i])
        }
        print(warnings())
        saveRDS(prob, file.)
    }
}

T.=T
pdf(file='TSVSlung.pdf')
plot(1:T., prob[ , 1], type='n', ylim=0:1, col=2,
     xlab='Steps', ylab='Inclusion Probability')
for(i in 1:P)
    if(prob[T., i]>0.5) {
        lines(1:T., prob[ , i], col=i)
        text(T.-i*10, prob[T.-i*10, i], dimnames(prob)[[2]][i],
             col=i, pos=2)
    }
dev.off()
sort(prob[T., ], TRUE)

