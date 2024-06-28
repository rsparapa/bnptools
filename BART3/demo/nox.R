
options(mc.cores=8)
library(BART3)
data(Boston) ## corrected data set

R=cor(Boston[ , -c(1:4, 7)])
R['nox', ]

x.train = Boston[ , -c(1:4, 7:8)]
post = mc.gbart(x.train, Boston$CMEDV, seed=99, sparse = TRUE)

print(cor(Boston$CMEDV, post$yhat.train.mean[-(1:10)])^2)

print(sort(post$varprob.mean, TRUE))

L=41
(x=seq(min(x.train[ , 'nox']), max(x.train[ , 'nox']), length.out=L))

print(system.time(pred <- FPD(post, x, S = 7))) ## assuming independence

plot(x, pred$yhat.test.mean, type='l', lwd=2,
     xlab='nox: Nitrogen Oxides air pollution',
     ylab='mdev: median home value (in thousands)',
     ylim=c(0, 50))
lines(x, pred$yhat.test.lower, lty=2, lwd=2)
lines(x, pred$yhat.test.upper, lty=2, lwd=2)
abline(h=c(0, 50), col='gray', lwd=2)
abline(h=c(mean(y)), col='blue', lwd=2)

## model similar to that presented in Harrison and Rubinfeld (1978)
fit=lm(log(CMEDV)~I(rm^2)+age+I(log(dis))+I(log(rad))+tax+ptratio+
           I((black-0.63)^2)+I(log(lstat))+crim+zn+indus+chas+
           I((nox-0.55)^2), data=Boston)
summary(fit)
lines(x, mean(Boston$CMEDV)*exp((x-0.55)^2*fit$coefficients[14]), lty=3, col='red', lwd=2)

dev.copy2pdf(file = 'nox.pdf')

##points(Boston$nox, Boston$CMEDV, col = 8)
