## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
##str(bmx)
##summary(bmx$RIDAGEEX)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

file.='growth1-fit1.rds'
if(!file.exists(file.)) {
    fit1 = mc.gbart(x.train, bmx$BMXHT, seed=21)
    saveRDS(fit1, file.)
} else { fit1=readRDS(file.) }
print(cor(bmx$BMXHT, fit1$yhat.train.mean)^2)

col.=c(4, 2) ## males=blue, females=red
plot(fit1$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
abline(b=1, a=0, col=8)

plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')

(N=length(bmx$BMXHT))
K=16 ##K=64
(age=seq(2, 18, length.out=K+1)[1:K])

x.test = fit1$x.train[1:(2*K), ]
x.test[ , 1]=rep(1:2, each=K)
x.test[ , 2]=age
print(x.test[ , c(1, 2)])

file.='growth1-pred1.rds'
if(!file.exists(file.)) {
    pred1 = FPD(fit1, x.test, c(1, 2))
    saveRDS(pred1, file.)
} else { pred1=readRDS(file.) }

pdf(file='growth1-pred1.pdf')
plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')
for(i in 1:2)
    lines(age, apply(pred1[ , (i-1)*K+1:K], 2, mean), col=col.[i], lwd=2)
dev.off()

age=c(age, 18)
K=length(age)
x.test = fit1$x.train[1:(2*K), ]
x.test[ , 1]=rep(1:2, each=K)
x.test[ , 2]=age

file.='growth1-pred2.rds'
if(!file.exists(file.)) {
    pred2 = FPD(fit1, x.test, 1:2, Subset=1:2)
    saveRDS(pred2, file.)
} else { pred2=readRDS(file.) }

pdf(file='growth1-pred2.pdf')
plot(bmx$RIDAGEEX, fit1$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')
for(i in 1:2) {
    lines(age, pred2$yhat.test.mean[(i-1)*K+1:K], col=col.[i], lwd=2)
    lines(age, pred2$yhat.test.lower[(i-1)*K+1:K], col=col.[i], lwd=2, lty=2)
    lines(age, pred2$yhat.test.upper[(i-1)*K+1:K], col=col.[i], lwd=2, lty=2)
}
dev.off()
