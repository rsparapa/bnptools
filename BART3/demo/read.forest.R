
library(BART3)
library(randomForest)

f = function(x)
    5+10*sin(pi*x[ , 1]*x[ , 2]) +  3*x[ , 3]^3

N = 2000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:P)
x.train[ , 4]= 1*(x.train[ , 4]>0)
y.train=(f(x.train)+sigma*rnorm(N))

H=100
x=seq(-2, 2, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

B=8
post1=mc.gbart(x.train, y.train, mc.cores=B, seed=12)

plot(post1$sigma[ , 1], type='l', ylim=c(0, max(c(post1$sigma))))
for(i in 2:B) lines(post1$sigma[ , i], col=i)
abline(h=sigma, lty=2)

post2=mc.gbart(x.train, y.train, mc.cores=B, seed=12, rfinit=TRUE)

for(i in 1:B) lines(post2$sigma[ , i], col=i, lty=2)

rf = randomForest(x.train, y.train, ntree=200, maxnodes=4)
getTree(rf, 1)

str(rf)

bMM = bartModelMatrix(x.train, numcut=100)

str(bMM$xinfo)

check = read.forest(rf, maxnodes=4, x.train=x.train)
write(check, 'check.txt')

