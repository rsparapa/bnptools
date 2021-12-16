
library(nftbart)

a = function(x) {
    x=rbind(x)
    exp(-2+4*(0.4*x[ , 1]+0.2*x[ , 2]-0.6*x[ , 2]*x[ , 3]))
}

b = function(x) {
    x=rbind(x)
    2.0-1.5*x[ , 4]+0.5*x[ , 5]+2*x[ , 5]*x[ , 6]
}

N=500
P.=c(6, 20)
odds.=list(seq(1, P.[1], 2), seq(1, P.[2], 2))
post=list()
y=list()
for(j in 1:2) {
    set.seed(42*j)
    P=P.[j]
    odds=odds.[[j]]
    x.train=matrix(runif(N*P), nrow=N, ncol=P)
    x.train[ , odds]=1*(x.train[ , odds]>0.5)
    dimnames(x.train)[[2]]=paste0('x', 1:P)
    y[[j]]=exp(rnorm(N, b(x.train), a(x.train)))
    C=exp(rnorm(N, b(x.train), a(x.train)))
    delta=1*(y[[j]]<C)
    times=delta*y[[j]]+(1-delta)*C
    set.seed(12*j)
    post[[j]]=nft(x.train, times, delta, K=0)
}

print(post[[1]]$lpml-post[[2]]$lpml)
print(post[[1]]$LPML-post[[2]]$LPML)

M=2000
r=list(0, 0)
for(j in 1:2) {
    z=log(y[[j]])
    for(i in 1:M) r[[j]][i]=cor(post[[j]]$z.train[i, ], z)
}
summary(r[[1]]^2)
summary(r[[2]]^2)
summary((r[[1]]^2)-(r[[2]]^2))
