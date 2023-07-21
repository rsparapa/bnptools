## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)
summary(bmx$RIDAGEEX)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

file.='growth1-fit1.rds'
if(!file.exists(file.)) {
    fit1 = mc.gbart(x.train, bmx$BMXHT, seed=21)
    saveRDS(fit1, file.)
} else { fit1=readRDS(file.) }
print(cor(bmx$BMXHT, fit1$yhat.train.mean)^2)

(N=length(bmx$BMXHT))
K=16 ##K=64
(age=seq(2, 18, length.out=K+1)[1:K])
x.test=cbind(rep(1:2, each=K), age)
print(x.test)

file.='growth2-fit0.rds'
if(!file.exists(file.)) {
    fit0 = mc.gbart(x.train[ , -(3:4)], bmx$BMXWT, x.test, seed=21)
    saveRDS(fit0, file.)
} else { fit0=readRDS(file.) }
print(cor(bmx$BMXWT, fit0$yhat.train.mean)^2)

P=ncol(fit1$x.train)
x.test = fit1$x.train[1:(2*K), ]
x.test[ , 1]=rep(1:2, each=K)
x.test[ , 2]=age
x.test[ , P]=fit0$yhat.test.mean
print(x.test[ , c(1, 2, P)])

pred0 = FPD(fit1, x.test, c(1, 2, P))
pred1 = FPDK(fit1, x.test, c(1, 2, P), kern.var=FALSE, mult.impute=5L)
##pred2 = FPDK(fit1, x.test, c(1, 2, P)) ## kern.var=TRUE default
pred3 = FPDK(fit1, x.test, c(1, 2, P), mult.impute=5L) ## kern.var=TRUE default

pdf(file='growth4-FPD.pdf')
plot(age, pred0$yhat.test.mean[1:K], col=4, lwd=2,
     type='l', xlab='Age', ylab='Height (cm)')
lines(age, pred0$yhat.test.mean[K+1:K], col=2, lwd=2)
for(i in 1:2) {
    points(age, pred1$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i))
    points(age, pred3$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i), pch='x')
}
legend('topleft', c('FPD', 'FPDK (K=5)', 'FPDK (K=30)'),
       pch=c(' ', 'x', 'o'), lty=c(1, 0, 0), lwd=c(2, 0, 0))
dev.off()

pdf(file='growth4-FPD-M.pdf')
plot(bmx$RIDAGEEX, bmx$BMXHT,
xlim=c(12, 17), ylim=c(160, 180),
     type='n', xlab='Age', ylab='Height (cm)',
sub='M only: 95% credible interval')
lines(age, pred0$yhat.test.lower[1:K], col=2, lwd=2)
lines(age, pred0$yhat.test.upper[1:K], col=4, lwd=2)
lines(age, pred1$yhat.test.lower[1:K], col=2, lwd=2, lty=2)
lines(age, pred1$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
lines(age, pred3$yhat.test.lower[1:K], col=2, lwd=2, lty=3)
lines(age, pred3$yhat.test.upper[1:K], col=4, lwd=2, lty=3)
legend('topleft', c('FPD', 'FPDK (w/o EV, K=5)', 'FPDK (w/ EV, K=5)'),
       lty=1:3, lwd=2)
dev.off()

source('check.R')

pred4 = SHAPK(fit1, x.test, c(1, 2, P), kern.var=FALSE, mult.impute=5L)
pred5 = SHAPK(fit1, x.test, c(1, 2, P)) 
pred6 = SHAPK(fit1, x.test, c(1, 2, P), mult.impute=5L)## kern.var=TRUE default
pred7 = SHAP(fit1, x.test, c(1, 2, P))
pred8 = SHAP(fit1, x.test, S=c(1, 2))

pdf(file='growth4-SHAPwt.pdf')
plot(age, pred0$yhat.test.mean[1:K], col=4, lwd=2,
     type='l', xlab='Age', ylab='Height (cm)')
lines(age, pred0$yhat.test.mean[K+1:K], col=2, lwd=2)
for(i in 1:2) {
    points(age, fit1$offset+pred7$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i), pch='x')
    points(age, fit1$offset+pred8$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i), pch='o')
}
legend('topleft', c('FPD (w/ wt)', 'mu+SHAP (w/ wt)', 'mu+SHAP (w/o wt)'),
       pch=c(' ', 'x', 'o'), lty=c(1, 0, 0), lwd=c(2, 0, 0))
dev.off()

pdf(file='growth4-SHAP.pdf')
plot(age, pred0$yhat.test.mean[1:K], col=4, lwd=2,
     type='l', xlab='Age', ylab='Height (cm)')
lines(age, pred0$yhat.test.mean[K+1:K], col=2, lwd=2)
for(i in 1:2) {
    points(age, fit1$offset+pred4$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i))
    ##points(age, fit1$offset+apply(pred9[ , (i-1)*K+1:K], 2, mean), col=2^(3-i), pch='x')
    points(age, fit1$offset+pred6$yhat.test.mean[(i-1)*K+1:K], col=2^(3-i), pch='x')
}
legend('topleft', c('FPD', 'mu+SHAP', 'mu+SHAPK (K=30)'),
       pch=c(' ', 'x', 'o'), lty=c(1, 0, 0), lwd=c(2, 0, 0))
dev.off()

pdf(file='growth4-SHAP-M.pdf')
plot(bmx$RIDAGEEX, bmx$BMXHT,
##xlim=c(12, 17), ylim=c(160, 180),
     type='n', xlab='Age', ylab='Height (cm)',
     sub='M only: 95% credible interval')
lines(age, pred0$yhat.test.lower[1:K], col=2, lwd=2)
lines(age, pred0$yhat.test.upper[1:K], col=4, lwd=2)
lines(age, fit1$offset+pred4$yhat.test.lower[1:K], col=2, lwd=2, lty=2)
lines(age, fit1$offset+pred4$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
lines(age, fit1$offset+pred5$yhat.test.lower[1:K], col=2, lwd=2, lty=3)
lines(age, fit1$offset+pred5$yhat.test.upper[1:K], col=4, lwd=2, lty=3)
legend('topleft', c('FPD', 'mu+SHAPK (w/ EV, K=5)', 'mu+SHAPK (w/ EV, K=30)'),
       lty=1:3, lwd=2)
dev.off()

