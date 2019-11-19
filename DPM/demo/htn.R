
library(DPM)

set.seed(12)

M = 100
N = 4*M
SD = 5
x = rt(3*M, 16)*SD+seq(110, 130, length.out=3*M)
y = rt(M,    4)*SD+seq(135, 150, length.out=M)

z = c(x, y)

##set.seed(21); fit <- mutau(z)
set.seed(21); fit <- mutau(z, b0.draw=0, b0=10000*SD^2)

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

hist(z, xlim=c(100, 160), breaks=25, freq=FALSE)
curve(density., add=TRUE, xlim=c(100, 160), n=201)

##dev.copy2pdf(file='dpmlio-galaxy.pdf')
