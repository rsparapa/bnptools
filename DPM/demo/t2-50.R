
library(DPM)

set.seed(12)
N <- 400
df=2
Y <- rt(N, df)
Z <- rt(N, df)
delta=(Y<Z)
y=delta*Y+(1-delta)*Z
table(delta)/N
A=-6
B=6
x=seq(A, B, length.out=1001)

set.seed(21); fit <- mutau(y, delta)
##tau.init=1/var(Y)
##set.seed(21); fit <- mutau(Y, a0=2, b0.b=tau.init, k0.b=tau.init*7.5)
##set.seed(21); fit <- mutau(Y, b0.draw=0)
plot(sapply(fit, function (x) x$hyper$alpha), type='l')
plot(sapply(fit, function (x) max(x$C)), type='l')
table(sapply(fit, function (x) max(x$C)))/length(fit)
plot(sapply(fit, function (x) x$hyper$k0), type='l')
plot(sapply(fit, function (x) log(x$hyper$k0)), type='l')
plot(sapply(fit, function (x) x$hyper$b0), type='l')
plot(sapply(fit, function (x) log(x$hyper$b0)), type='l')

pdf. <- function(x)
    apply(sapply(fit, function(a) {
        b <- 0
        for(k in 1:max(a$C))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, mean)

pdf.025 <- function(x)
    apply(sapply(fit, function(a) {
        b <- 0
        for(k in 1:max(a$C))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=0.025)

pdf.975 <- function(x)
    apply(sapply(fit, function(a) {
        b <- 0
        for(k in 1:max(a$C))
            b <- b+a$states[k]*dnorm(x, mean=a$phi[k, 1],
                                     sd=1/sqrt(a$phi[k, 2]))
        return(b/N)}), 1, quantile, probs=0.975)

plot(x, dt(x, df), col='blue', ylab='f(x)', type='l', lwd=2)
##lines(x, dt(x, df), col='red', lwd=2)
lines(x, pdf.025(x), col='green', lwd=2, lty=2)
lines(x, pdf.(x), col='green', lwd=2)
lines(x, pdf.975(x), col='green', lwd=2, lty=2)
rug(y)
abline(v=0, h=0, col='gray')
legend('topright',
       legend=c('t(2)', 'DPM LIO Neal 7'),
       lty=1, col=c('blue', 'green'), lwd=2)
dev.copy2pdf(file='t2-50-good.pdf')

