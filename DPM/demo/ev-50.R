
library(DPM)
library(extraDistr)

set.seed(12)
N <- 500
scale=0.5
mu=5
Y <- rgev(N, mu, sigma=scale)
Z <- rgev(N, mu, sigma=scale)
delta=(Y<Z)
y=delta*Y+(1-delta)*Z
table(delta)/N
A=-6*scale
x=mu+seq(A, A+12*scale, length.out=1001)

set.seed(21); fit7 <- mutau(y, delta, neal.m=0)
set.seed(21); fit8 <- mutau(y, delta, neal.m=1)

plot(sapply(fit8, function (x) x$hyper$alpha), type='l')
plot(sapply(fit8, function (x) max(x$C)), type='l')
table(sapply(fit8, function (x) max(x$C)))/length(fit8)
plot(sapply(fit8, function (x) x$hyper$k0), type='l')
plot(sapply(fit8, function (x) log(x$hyper$k0)), type='l')
plot(sapply(fit8, function (x) x$hyper$b0), type='l')
plot(sapply(fit8, function (x) log(x$hyper$b0)), type='l')

pdf7CI <- function(x)
    apply(sapply(fit7, function(a) {
        b <- 0
        for(k in 1:max(a$C))
        ##for(k in seq_along(a$states))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

pdf7=pdf7CI(x)

pdf8CI <- function(x)
    apply(sapply(fit8, function(a) {
        b <- 0
        for(k in 1:max(a$C))
        ##for(k in seq_along(a$states))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=c(0.025, 0.5, 0.975))

pdf8=pdf8CI(x)

plot(x, dgev(x, mu, scale), col='blue', ylab='f(x)', type='l', lwd=2,
     ylim=c(0, 1.2*dgev(mu, mu, scale)),
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
       legend=c(paste0('EV(', mu, ', ', scale, ')'),
                'DPM Neal7', 'DPM Neal8'),
       lty=1, col=c('blue', 'red', 'green'), lwd=2)
##dev.copy2pdf(file='t1-50-NS.pdf')



## fit=mutau(Y)

## fit <- survreg(Surv(time=exp(y), event=delta)~1,
##                data=data.frame(y=y, delta=delta))
## plot(x, dgev(x, mu, scale), col='blue', ylab='f(x)', type='l')
## lines(x, dgev(x, fit$coef-fit$scale*0.5772, fit$scale), col='red')
## lines(x, dgev(x, fit$coef, fit$scale), col='black')
## summary(fit)

## pgev(mu, mu, scale)
## pgev(mu-scale*0.5772, mu, scale)

est.q <- function(p, fit)
    mean(sapply(fit, function(a) {
        b <- 0
        for(k in 1:max(a$C))
            b <- b+a$states[k]*qnorm(p, a$phi[k, 1], a$phi[k, 2]^(-0.5))
        return(b/N)}))
est.q(0.05, fit8)
qgev(0.05, mu, scale)
est.q(0.5, fit8)
qgev(0.5, mu, scale)
est.q(0.95, fit8)
qgev(0.95, mu, scale)

