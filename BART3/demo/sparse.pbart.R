
library(BART3)

B <- 8 ##getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-0.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1  #y = f(x) + sigma*z where z~N(0, 1)
P = 100    #number of covariates
##thin <- c(10, 10, 10)
n <- c(250, 500, 1000, 2000, 4000, 8000)
M=length(n)
N <- max(n)
set.seed(12)
x.train=matrix(runif(N*P), N, P)
y.train=(rnorm(N, f(x.train), sigma)>0)*1

file.='sparse-pbart.rds'
## if(file.exists(file.)) {
##     post=readRDS(file.)
## } else
{
    post <- list()

    for(i in 1:M) {
        N=n[i]
        post[[i]] = mc.gbart(x.train[1:N, ], y.train[1:N],
                             mc.cores=B, type='pbart',
                             sparse=TRUE, seed=99)
    }
    saveRDS(post, file.)
}

par(mfrow=c(3, 1))

for(i in 1:M) {
    N <- n[i]
    
    if(N %in% c(250, 1000, 4000)) {
        plot(post[[i]]$varprob.mean, col=c(rep(2, 5), rep(1, P-5)),
             sub=paste0('N:', N, ', P:', P), ##, ', thin:', thin[i]),
             ylab='Selection Probability', ylim=c(0, 0.3),
             pch=1+45*(post[[i]]$varprob.mean <= 1/P))
        lines(c(0, 100), c(1/P, 1/P))

        table(1+45*(post[[i]]$varprob.mean <= 1/P))
    }
}

par(mfrow=c(1, 1))

if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'sparse-pbart.pdf', sep='/'))
dev.copy2pdf(file='sparse-pbart.pdf')

for(i in M:1) {    
    N <- n[i]
    if(i==M) vimp1=c(N, sum(post[[i]]$varprob.mean > 1/P))
    else vimp1=rbind(vimp1, c(N, sum(post[[i]]$varprob.mean > 1/P)))
    if(i==M) {
        vimp3=1*(post[[i]]$varprob.mean > 1/P)
        vimp2=c(N, sum(vimp3))
    }
    else {
        vimp3=vimp3*(post[[i]]$varprob.mean > 1/P)
        vimp2=rbind(vimp2, c(N, sum(vimp3)))
    }
}
print(vimp1)
print(vimp2)
