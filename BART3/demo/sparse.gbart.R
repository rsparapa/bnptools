
B = 8 ## number of threads
options(mc.cores = B)
library(BART3)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) ## only the first 5 matter
    10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-0.5)^2 + 10*x[ , 4] + 5*x[ , 5]

P = 125                   ## number of covariates: signal+noise
C = 60                    ## number of continuous: noise
C. = C+5                  ## number of continuous: signal+noise
D = P-C.                  ## number of dichotomous noise
nskip = 1000              ## burn-in discarded
ndpost = 1000             ## total length of combined chain
ndpost. = floor(ndpost/B) ## length of each chain by thread

par(mfrow=c(2, 2))

post <- list()

for(i in 1:2) {
    n <- c(100, 500, 2500)[i]
    set.seed(34)
    x.train=matrix(c(runif(n*C.), rbinom(n*D, 1, 0.5)), n, P)
    y.train=rnorm(n, f(x.train))

    for(j in c(TRUE, FALSE)) {
        h <- (i-1)*2+j+1
        post[[h]] = mc.gbart(x.train, y.train, sparse=j, 
                             seed=99, ndpost=ndpost, nskip=nskip)
        
        if(j) {
            varprob <- post[[h]]$varprob.mean
            varprob. <- post[[h]]$varprob[1:ndpost., 1]
        } else {
            varprob <- post[[h]]$varcount.mean/sum(post[[h]]$varcount.mean)
            varprob. <- apply(post[[h]]$varcount, 1, sum)
            varprob. <- (post[[h]]$varcount[ , 1]/varprob.)[1:ndpost.]
        }

        col = c(rep(2, 5), rep(1, C), rep(4, D))
        plot(varprob, col=col, pch = '.',
             main=paste0('N:', n, ', P:', P, ', DART:', j),
             ylab='Selection Probability', ylim=0:1)
        if(j) {
            lines(P*(1:ndpost.)/ndpost., varprob., col = 8, lty = 3)
            text(P, varprob[1], expression(s[1]))
        }
        points(varprob, col=col, pch = c('.', 'o')[1+(varprob>1/P)])
        if(i == j) 
            legend('topleft', 
               c('Signal: uniform', 'Noise: uniform', 'Noise: Bernoulli'),
               col = c(2, 1, 4), pch = 'o')
    }
}

par(mfrow=c(1, 1))
dev.copy2pdf(file = 'sparse-gbart.pdf')
