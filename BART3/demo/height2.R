## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)
##summary(bmx$RIDAGEEX)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

## fit1 = mc.gbart(x.train, bmx$BMXHT, seed=21)
## print(sort(fit1$varprob.mean, TRUE))
## print(cor(bmx$BMXHT, fit1$yhat.train.mean)^2)

## col.=c(4, 2) ## males=blue, females=red
## plot(fit1$yhat.train.mean, bmx$BMXHT, asp=1,
##      pch='.', col=col.[bmx$RIAGENDR],
##      xlab='Predicted Height (cm)',
##      ylab='Observed Height (cm)')
## abline(b=1, a=0, col=8)

fit2 = mc.gbart(x.train, bmx$BMXHT, seed=21, sparse=TRUE)
print(sort(fit2$varprob.mean, TRUE))
print(cor(bmx$BMXHT, fit2$yhat.train.mean)^2)

col.=c(4, 2) ## males=blue, females=red
plot(fit2$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
abline(b=1, a=0, col=8)

plot(bmx$RIDAGEEX, fit2$yhat.train.mean,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')

(N=length(bmx$BMXHT))
K=64
(age=seq(2, 18, length.out=K+1)[1:K])

x.test = x.train
for(k in 2:K) x.test=rbind(x.test, x.train)
x.test$RIDAGEEX=rep(age, each=N)
str(x.test)

pred2 = predict(fit2, x.test)

marg2 = 0
marg2.025 =0
marg2.975 =0
for(k in 1:K) {
    marg2[k] = mean(apply(pred2[ , (k-1)*N+1:N], 1, mean))
    marg2.025[k] = quantile(apply(pred2[ , (k-1)*N+1:N], 1, mean), probs=0.025)
    marg2.975[k] = quantile(apply(pred2[ , (k-1)*N+1:N], 1, mean), probs=0.975)
}

lines(age, marg2, lwd=2)
lines(age, marg2.025, lty=2)
lines(age, marg2.975, lty=2)

x.train.=bmx[ , 3:4]
print(cor(bmx$BMXWT, x.train.[ , 2])^2)
fit0 = mc.gbart(x.train., bmx$BMXWT, seed=20)
print(cor(bmx$BMXWT, fit0$yhat.train.mean)^2)

x.test. = x.train.
for(k in 2:K) x.test.=rbind(x.test., x.train.)
x.test.$RIDAGEEX=rep(age, each=N)
x.test.=rbind(x.test., x.test.)
x.test.$RIAGENDR=rep(1:2, each=N*K)
str(x.test.)

pred0 = predict(fit0, x.test.)

marg0 = 0
for(k in 1:(2*K)) marg0[k] = mean(apply(pred0[ , (k-1)*N+1:N], 1, mean))

x.test0 = x.train
for(k in 2:K) x.test0 = rbind(x.test0, x.train)
x.test0$RIDAGEEX=rep(age, each=N)
x.test1=x.test0
x.test1$RIAGENDR=1
x.test1$BMXWT=rep(marg0[1:K], each=N)
x.test2=x.test0
x.test2$RIAGENDR=2
x.test2$BMXWT=rep(marg0[K+1:K], each=N)
x.test0=rbind(x.test1, x.test2)
str(x.test0)

pred2. = predict(fit2, x.test0)

marg2. = 0
for(k in 1:(2*K)) marg2.[k] = mean(apply(pred2.[ , (k-1)*N+1:N], 1, mean))
lines(age, marg2.[1:K], lwd=2, col=4)
lines(age, marg2.[K+1:K], lwd=2, col=2)
##dev.copy2pdf(file='bart-growth.pdf')

                     
