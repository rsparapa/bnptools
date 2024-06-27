
options(mc.cores=8)
library(BART3)
data(Boston) ## corrected data set

(B <- getOption('mc.cores', 1))
figures = getOption('figures', default='NONE')

y = Boston$CMEDV ## median value
x.train = Boston[ , -c(1:4, 7:8)]
(N=length(y))   ## total sample size
post = mc.gbart(x.train, y, mc.cores=B, seed=99)

R=cor(cbind(y, x.train))
R['nox', ]

L=41
(x=seq(min(x.train[ , 'nox']), max(x.train[ , 'nox']), length.out=L))

print(system.time(pred <- FPD(post, x, S = 7))) ## assuming independence

plot(x, pred$yhat.test.mean, type='l', lwd=2,
     xlab='nox: Nitrogen Oxides air pollution',
     ylab='mdev: median home value (in thousands)',
     ylim=c(0, 50))
lines(x, pred$yhat.test.lower, lty=2, lwd=2)
lines(x, pred$yhat.test.upper, lty=2, lwd=2)
abline(h=c(0, 50), col='gray')

## model similar to that presented in Harrison and Rubinfeld (1978)
fit=lm(log(y)~I(rm^2)+age+I(log(dis))+I(log(rad))+tax+ptratio+
           I((black-0.63)^2)+I(log(lstat))+crim+zn+indus+chas+
           I((nox-0.55)^2), data=Boston)
summary(fit)
lines(x, mean(y)*exp((x-0.55)^2*fit$coefficients[14]), lty=3, col='red', lwd=2)

if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'nox.pdf', sep='/'))

if(figures!='NONE')
    dev.copy2eps(file=paste(figures, 'nox.eps', sep='/'))
