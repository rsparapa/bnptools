## variable selection example
options(mc.cores=8)
library(BART3)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)

x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

fit1 = mc.gbart(x.train, bmx$BMXHT, seed=21)
print(sort(fit1$varprob.mean, TRUE))
print(cor(bmx$BMXHT, fit1$yhat.train.mean)^2)

col.=c(4, 2) ## males=blue, females=red
plot(fit1$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
abline(b=1, a=0, col=8)

fit2 = mc.gbart(x.train, bmx$BMXHT, seed=21, sparse=TRUE)
print(sort(fit2$varprob.mean, TRUE))
print(cor(bmx$BMXHT, fit2$yhat.train.mean)^2)

col.=c(4, 2) ## males=blue, females=red
plot(fit2$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
abline(b=1, a=0, col=8)

(N=length(bmx$BMXHT))

## create knockoffs of weight
set.seed(20)
K=100
rho=0.5
##rho=0.85
Z=as.data.frame(matrix(0, nrow=N, ncol=K))
(mu.=mean(bmx$BMXWT))
(sd.=sd(bmx$BMXWT))
for(k in 1:K) Z[ , k]=rnorm(N, rho*(bmx$BMXWT-mu.), sqrt(1-(rho^2))*sd.)
print(cor(bmx$BMXWT, bmx$BMXHT))
print(cor(bmx$BMXWT, Z[ , 1:3]))
print(cor(bmx$BMXHT, Z[ , 1:3]))

fit3 = mc.gbart(cbind(x.train, Z), bmx$BMXHT, seed=22, sparse=TRUE)
(P=length(fit3$varprob.mean))
print(sort(fit3$varprob.mean[fit3$varprob.mean>(2/P)], TRUE))
print(cor(bmx$BMXHT, fit3$yhat.train.mean)^2)
    
