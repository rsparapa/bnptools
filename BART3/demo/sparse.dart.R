
library(BART3)

## the DART default with respect to variable selection probability
## weighs a category as much as a single continuous variable
## if it is coded as a factor; otherwise, each category is counted

f = function(x) #only the first 5 matter
    10*(sin(pi*x[ , 1]*x[ , 4]) + 2*(x[ , 5]-.5)^2+x[ , 2]+0.5*x[ , 3]-0.5)
sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
p = 100      #number of covariates
k=p-3
C = 8
n = 5000
set.seed(12)
## let 1, 2, 3 be categories
x.train=cbind(t(rmultinom(n, 1, rep(1, 3)/3)), matrix(runif(n*k), n, k))
dimnames(x.train)[[2]]=paste0('x', 1:p)
y.train = f(x.train)+sigma*rnorm(n)

par(mfrow=c(2, 1))
post=list()
for(i in 1:2) {
    if(i==1) X=x.train
    else {
        x.123 = x.train[ , 1:3] %*% c(1:3)
        X=data.frame(x1=factor(x.123))
        X=cbind(X, x.train[ , -c(1:3)])
    }

    post[[i]] = mc.gbart(X, y.train, mc.cores=C, sparse=TRUE, seed=99)

    plot(post[[i]]$varprob.mean, col=c(rep(2, 5), rep(1, p-5)),
         pch=1+44*(post[[i]]$varprob.mean<=1/p),
         ylab='Selection Probability', ylim=c(0, 1))
    lines(c(0, 100), c(1/p, 1/p))
    abline(v=3.5)
}
par(mfrow=c(1, 1))

check1=bartModelMatrix(x.train, numcut=100)
post[[1]]$grp
post[[1]]$rho
check2=bartModelMatrix(X, numcut=100)
post[[2]]$grp
post[[2]]$rho

plot(post[[1]]$varprob.mean, col=2, xlim=c(1, 5),
     ylab='Selection Probability', ylim=c(0.05, 0.45))
points(post[[2]]$varprob.mean, col=4, xlim=c(1, 5))
lines(c(0, 100), c(1/p, 1/p))
abline(v=3.5)
legend('topleft', c('matrix', 'data.frame/factor'),
       col=c(2, 4), pch=1)

