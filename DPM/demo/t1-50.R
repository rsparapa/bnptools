
library(DPM)

set.seed(12)
N <- 200
df=1
Y <- rt(N, df)
Z <- rt(N, df)
delta=(Y<Z)
y=delta*Y+(1-delta)*Z
table(delta)/N
A=-6
B=6
x=seq(A, B, length.out=1001)

set.seed(21); fit7 <- mutau(y, delta)
set.seed(21); fit8 <- mutau(y, delta, neal.m=1)

##tau.init=1/var(Y)
##set.seed(21); fit <- mutau(Y, a0=2, b0.b=tau.init, k0.b=tau.init*7.5)
##set.seed(21); fit <- mutau(Y, b0.draw=0)
plot(sapply(fit8, function (x) x$hyper$alpha), type='l')
plot(sapply(fit8, function (x) max(x$C)), type='l')
table(sapply(fit8, function (x) max(x$C)))/length(fit8)
plot(sapply(fit8, function (x) x$hyper$k0), type='l')
plot(sapply(fit8, function (x) log(x$hyper$k0)), type='l')
plot(sapply(fit8, function (x) x$hyper$b0), type='l')
plot(sapply(fit8, function (x) log(x$hyper$b0)), type='l')

## pdf. <- function(x)
##     apply(sapply(fit, function(a) {
##         b <- 0
##         for(k in seq_along(a$states))
##             b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
##                                      sd=1/sqrt(a$phi[k, 2]))
##         return(b/N)}), 1, mean)

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
     ylim=c(0, dt(0, df)),
     main=paste0('N=', N, ', 50% censoring, LIO non-standardized settings'))
##lines(x, dt(x, df), col='green', lwd=2)
lines(x, pdf7[1, ], col='red', lwd=2, lty=2)
lines(x, pdf7[2, ], col='red', lwd=2)
lines(x, pdf7[3, ], col='red', lwd=2, lty=2)
lines(x, pdf8[1, ], col='green', lwd=2, lty=2)
lines(x, pdf8[2, ], col='green', lwd=2)
lines(x, pdf8[3, ], col='green', lwd=2, lty=2)
rug(y)
abline(v=0, h=0, col='gray', lwd=2)
legend('topright',
       legend=c(paste0('t(', df, ')'), 'DPM Neal7', 'DPM Neal8'),
       lty=1, col=c('blue', 'red', 'green'), lwd=2)
##dev.copy2pdf(file='t1-50-NS.pdf')

cdf7CI <- function(x)
    apply(sapply(fit7, function(a) {
        b <- 0
        for(k in seq_along(a$states))
            b <- b+a$states[k]*pnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]),
                                     lower.tail=FALSE)
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

cdf7=cdf7CI(x)

cdf8CI <- function(x)
    apply(sapply(fit8, function(a) {
        b <- 0
        for(k in seq_along(a$states))
            b <- b+a$states[k]*pnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]),
                                     lower.tail=FALSE)
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

cdf8=cdf8CI(x)

km=survfit(Surv(time=exp(y), event=delta)~1)

X=exp(x)
plot(X, pt(x, df, lower.tail=FALSE), col='blue', ylab='S(t)', type='l', lwd=2,
     ylim=c(0, 1), xlim=c(0, 90), xlab='t=exp(y)',
     main=paste0('N=', N, ', 50% censoring, LIO non-standardized settings'))
lines(km, type='s', lwd=2)
##lines(x, dt(x, df), col='green', lwd=2)
lines(X, cdf7[1, ], col='red', lwd=2, lty=2)
lines(X, cdf7[2, ], col='red', lwd=2)
lines(X, cdf7[3, ], col='red', lwd=2, lty=2)
lines(X, cdf8[1, ], col='green', lwd=2, lty=2)
lines(X, cdf8[2, ], col='green', lwd=2)
lines(X, cdf8[3, ], col='green', lwd=2, lty=2)
abline(v=0, h=0:1, col='gray', lwd=2)
rug(y)
legend(65.25, 0.9975,
       legend=c(paste0('t(', df, ')'), 'Kaplan-Meier', 'DPM Neal7',
                'DPM Neal8'),
       lty=1, col=c('blue', 'black', 'red', 'green'), lwd=2)
##dev.copy2pdf(file='t1-50-NS-cdf.pdf')
