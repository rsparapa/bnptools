
library(DPM)
data(galaxy)

Y=c(galaxy$speed/1000)
N <- length(Y)

set.seed(21); fit <- mutau(Y)

plot(fit$Cmax, type='l')
table(fit$Cmax)/fit$mcmc$keep

plot(fit$alpha, type='l')
abline(h=1, col='gray')

density. <- function(x) {
    a = 0
    for(i in 1:fit$mcmc$keep) 
        for(j in 1:fit$Cmax[i]) {
            b = dnorm(x, fit$phi[i, j, 1], 1/sqrt(fit$phi[i, j, 2]))/N
            a = a+fit$states[i, j]*b/fit$mcmc$keep
        }
    return(a)
}

hist(Y, xlim=c(7, 37), breaks=25, freq=FALSE)
curve(density., add=TRUE, xlim=c(7, 37), n=201)

set.seed(21); fit2 <- NoGa(Y)

plot(fit2$Cmax, type='l')
table(fit2$Cmax)/fit2$mcmc$keep

plot(fit2$alpha, type='l')
abline(h=1, col='gray')

density. <- function(x) {
    a = 0
    for(i in 1:fit2$mcmc$keep) 
        for(j in 1:fit2$Cmax[i]) {
            b = dnorm(x, fit2$phi[i, j, 1], 1/sqrt(fit2$phi[i, j, 2]))/N
            a = a+fit2$states[i, j]*b/fit2$mcmc$keep
        }
    return(a)
}

hist(Y, xlim=c(7, 37), breaks=25, freq=FALSE)
curve(density., add=TRUE, xlim=c(7, 37), n=201)

##dev.copy2pdf(file='dpmlio-galaxy.pdf')
