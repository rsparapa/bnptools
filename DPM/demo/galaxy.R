
library(DPM)
data(galaxy)

Y=c(galaxy$speed/1000)
N <- length(Y)

##set.seed(21); fit <- mutau(Y)
set.seed(21); fit <- mutau(Y, neal.m=1)
##set.seed(21); fit <- mutau(Y, keep=3, burn=0, thin=1, neal.m=1)
        
##hist(sapply(fit, function (x) max(x$C))+1, freq=FALSE)

plot(sapply(fit, function (x) x$hyper$alpha), type='l')
abline(h=1, col='gray')

plot(sapply(fit, function (x) max(x$C)), type='l')

table(sapply(fit, function (x) max(x$C)))/length(fit)

density. <- function(x)
    apply(sapply(fit, function(a) {
        b <- 0
        for(k in seq_along(a$states))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, mean)

hist(Y, xlim=c(7, 37), breaks=25, freq=FALSE)
curve(density., add=TRUE, xlim=c(7, 37), n=201)

##dev.copy2pdf(file='dpmlio-galaxy.pdf')
