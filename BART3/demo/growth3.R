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

pdf(file='growth1-zep1.pdf')
col.=c(4, 2) ## males=blue, females=red
plot(fit1$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
abline(b=1, a=0, col=8)
dev.off()

pdf(file='growth1-fit1.pdf')
plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')
dev.off()

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

pred1 = FPDK(fit1, x.test, c(1, 2, P), kern.var=FALSE)
pred2 = FPDK(fit1, x.test, c(1, 2, P)) ## kern.var=TRUE default

pdf(file='growth3-FPDK.pdf')

par(mfrow=c(1, 2))
plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR], type='n',
     ylab='Height (cm)', xlab='Age (yr)', sub='kern.var=FALSE')
for(i in 1:2) {
    lines(age, pred1$yhat.test.mean[(i-1)*K+1:K], col=col.[i])
    lines(age, pred1$yhat.test.lower[(i-1)*K+1:K],
          lty=2, col=col.[i], lwd=2)
    lines(age, pred1$yhat.test.upper[(i-1)*K+1:K],
          lty=2, col=col.[i], lwd=2)
}

plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR], type='n',
     ylab='', xlab='Age (yr)', sub='kern.var=TRUE')
for(i in 1:2) {
    lines(age, pred2$yhat.test.mean[(i-1)*K+1:K], col=col.[i])
    lines(age, pred2$yhat.test.lower[(i-1)*K+1:K],
          lty=2, col=col.[i], lwd=2)
    lines(age, pred2$yhat.test.upper[(i-1)*K+1:K],
          lty=2, col=col.[i], lwd=2)
}
par(mfrow=c(1, 1))
dev.off()
