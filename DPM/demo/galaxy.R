
library(DPM)
data(galaxy)

Y=c(galaxy$speed/1000)
N <- length(Y)

set.seed(21); fit <- mutau(Y)

Cmax = apply(fit$C, 1, max)
plot(Cmax, type='l')
table(Cmax)/fit$mcmc$keep

plot(fit$alpha, type='l')
abline(h=1, col='gray')

density. <- function(x) {
    a = 0
    for(i in 1:fit$mcmc$keep) 
        for(j in 1:Cmax[i]) {
            b = dnorm(x, fit$phi[i, j, 1], 1/sqrt(fit$phi[i, j, 2]))/N
            a = a+fit$states[i, j]*b/fit$mcmc$keep
        }
    return(a)
}

hist(Y, xlim=c(7, 37), breaks=25, freq=FALSE)
curve(density., add=TRUE, xlim=c(7, 37), n=201)

##dev.copy2pdf(file='dpmlio-galaxy.pdf')
