
library(BART3)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) ##only the first 5 matter
    10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-0.5)^2 + 10*x[ , 4] + 5*x[ , 5]

k = 100      ## number of covariates
c = 48       ## number of continuous noise variables
c. = c+5     ## number of continuous variables
d = k-c.     ## number of dichotomous noise variables
ndpost = 1000
nskip = 1000
C = 8

par(mfrow=c(2, 2))

post <- as.list(1:6)

for(i in 1:2) {
    n <- c(100, 500, 2500)[i]
    set.seed(34)
    x.train=matrix(c(runif(n*c.), rbinom(n*d, 1, 0.5)), n, k)
    y.train=rnorm(n, f(x.train))

    for(j in c(TRUE, FALSE)) {
        h <- (i-1)*2+j+1
        post[[h]] = mc.gbart(x.train, y.train, mc.cores=C, sparse=j, 
                             seed=99, ndpost=ndpost, nskip=nskip)
        
        if(j) {
            varprob <- post[[h]]$varprob.mean
            varprob. <- post[[h]]$varprob[ , 1]
        } else {
            varprob <- post[[h]]$varcount.mean/sum(post[[h]]$varcount.mean)
            varprob. <- apply(post[[h]]$varcount, 1, sum)
            varprob. <- post[[h]]$varcount[ , 1]/varprob.
        }

        col = c(rep(2, 5), rep(1, c), rep(4, d))
        plot(varprob, col=col, pch = '.',
             main=paste0('N:', n, ', P:', k, ', DART:', j),
             ylab='Selection Probability', ylim=0:1)
        abline(h = 1/k)
        lines(k*(1:ndpost)/ndpost, varprob., col = 8, lty = 3)
        points(varprob, col=col, pch = c('.', 'o')[1+(varprob>1/k)])
        if(i == j) 
            legend('topleft', 
               c('Signal: uniform', 'Noise: uniform', 'Noise: Bernoulli'),
               col = c(2, 1, 4), pch = 'o')
    }
}

par(mfrow=c(1, 1))

