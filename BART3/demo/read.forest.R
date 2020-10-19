
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

rf = randomForest(x.train, y.train, forest=TRUE, ntree=200, maxnodes=4)
getTree(rf, 1)

str(rf)

bMM = bartModelMatrix(x.train, numcut=100)

str(bMM$xinfo)

check = read.forest(rf, maxnodes=4, x.train=x.train)
write(check, 'check.txt')

## h=!is.na(c(check[[1]][ , , 3]))
## cor(c(check[[1]][ , , 3])[h], c(check[[1]][ , , 5])[h])
## summary(c(check[[1]][ , , 3])[h])
## summary(c(check[[1]][ , , 5])[h])
## summary(abs(c(check[[1]][ , , 3])[h]-c(check[[1]][ , , 5])[h]))

## nodes = integer(200)
## for(i in 1:200) nodes[i]=nrow(getTree(rf, i))
## print(table(nodes))
## for(i in which(nodes<7)) print(dimnames(getTree(rf, i))[[1]])
