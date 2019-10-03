
library(DPM)

set.seed(12)
N <- 100
df=1
Y <- rt(N, df)
A=-6
B=6
x=seq(A, B, length.out=1001)

dpm.lio <- normwrapper(y=Y, ngrid=1001, grid=x)

set.seed(21); fit8 <- mutau(Y, neal.m=1)
set.seed(21); fit7 <- mutau(Y)
##tau.init=1/var(Y)
## plot(sapply(fit, function (x) x$hyper$alpha), type='l')
## plot(sapply(fit, function (x) max(x$C)), type='l')
## table(sapply(fit, function (x) max(x$C)))/length(fit)
## plot(sapply(fit, function (x) x$hyper$k0), type='l')
## plot(sapply(fit, function (x) log(x$hyper$k0)), type='l')
## plot(sapply(fit, function (x) x$hyper$b0), type='l')
## plot(sapply(fit, function (x) log(x$hyper$b0)), type='l')

pdf7CI <- function(x)
    apply(sapply(fit7, function(a) {
        b <- 0
        for(k in seq_along(a$states))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

pdf7=pdf7CI(x)

pdf8CI <- function(x)
    apply(sapply(fit8, function(a) {
        b <- 0
        for(k in seq_along(a$states))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

pdf8=pdf8CI(x)

plot(x, dt(x, df), col='blue', ylab='f(x)', type='l', lwd=2,
     ylim=c(0, dt(0, 2)),
     main=paste0('N=', N, ', no censoring, LIO non-standardized settings'))
lines(dpm.lio$dens~dpm.lio$grid1, type='l', col='black', lwd=2)
lines(x, pdf7[1, ], col='red', lwd=2, lty=2)
lines(x, pdf7[2, ], col='red', lwd=2)
lines(x, pdf7[3, ], col='red', lwd=2, lty=2)
lines(x, pdf8[1, ], col='green', lwd=2, lty=2)
lines(x, pdf8[2, ], col='green', lwd=2)
lines(x, pdf8[3, ], col='green', lwd=2, lty=2)
rug(Y)
abline(v=0, h=0, col='gray', lwd=2)
legend('topright',
       legend=c(paste0('t(', df, ')'),
                'DPpackage Neal8', 'DPM Neal7', 'DPM Neal8'),
       lty=1, col=c('blue', 'black', 'red', 'green'), lwd=2)
##dev.copy2pdf(file='t1-0-NS.pdf')

