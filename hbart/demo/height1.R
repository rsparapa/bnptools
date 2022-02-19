## marginal effects example
options(mc.cores=8)
library(hbart)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

set.seed(21)
X.train=hbartModelMatrix(x.train)
fit2 = hbart(X.train, bmx$BMXHT)
print(sort(fit2$mu.varprob, TRUE))
print(cor(bmx$BMXHT, fit2$pred$f.train.mean)^2)

col.=c(4, 2) ## males=blue, females=red
## plot(fit2$pred$f.train.mean, bmx$BMXHT, asp=1,
##      pch='.', col=col.[bmx$RIAGENDR],
##      xlab='Predicted Height (cm)',
##      ylab='Observed Height (cm)')
## abline(b=1, a=0, col=8)

(N=length(bmx$BMXHT))
K=64
(age=seq(2, 18, length.out=K+1)[1:K])

x.train.=bmx[ , 3:4]
print(cor(bmx$BMXWT, x.train.[ , 2])^2)
set.seed(20)
X.train. = hbartModelMatrix(x.train.)
fit0 = hbart(X.train., bmx$BMXWT)
print(cor(bmx$BMXWT, fit0$pred$f.train.mean)^2)

x.test. = x.train.
for(k in 2:K) x.test.=rbind(x.test., x.train.)
x.test.$RIDAGEEX=rep(age, each=N)
x.test.=rbind(x.test., x.test.)
x.test.$RIAGENDR=rep(1:2, each=N*K)
str(x.test.)

X.test. = hbartModelMatrix(x.test.)
pred0 = predict(fit0, X.test.)

pred0.= pred0$f.test
(M = nrow(pred0.))
(L = ncol(pred0.)/(N*K))
marg0 = 0
for(l in 1:L) { 
    h=(l-1)*K
    for(k in 1:K) {
        j=h+k
        marg0[j] = mean(apply(pred0.[ , (j-1)*N+1:N], 1, mean))
    }
}

marg0=c(sort(marg0[1:K]), sort(marg0[K+1:K]))

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

X.test0 =hbartModelMatrix(x.test0)
pred2 = predict(fit2, X.test0)
pred2.= pred2$f.test
marg2. = 0
marg2.025 = 0
marg2.975 = 0
for(k in 1:(2*K)) {
    fpd = apply(pred2.[ , (k-1)*N+1:N], 1, mean)
    marg2.[k] = mean(fpd)
    marg2.025[k] = quantile(fpd, probs=0.025)
    marg2.975[k] = quantile(fpd, probs=0.975)
}

pdf('nosort-growth.pdf')
plot(bmx$RIDAGEEX, bmx$BMXHT,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age, marg2.[1:K], lwd=2, col=4)
lines(age, marg2.025[1:K], lty=2, col=4)
lines(age, marg2.975[1:K], lty=2, col=4)
lines(age, marg2.[K+1:K], lwd=2, col=2)
lines(age, marg2.025[K+1:K], lty=2, col=2)
lines(age, marg2.975[K+1:K], lty=2, col=2)
dev.off()
##dev.copy2pdf(file='nosort-growth.pdf')

marg2.=c(sort(marg2.[1:K]), sort(marg2.[K+1:K]))
marg2.025=c(sort(marg2.025[1:K]), sort(marg2.025[K+1:K]))
marg2.975=c(sort(marg2.975[1:K]), sort(marg2.975[K+1:K]))

pdf('sort-growth.pdf')
plot(bmx$RIDAGEEX, bmx$BMXHT,
     pch='.', col=col.[bmx$RIAGENDR],
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age, marg2.[1:K], lwd=2, col=4)
lines(age, marg2.025[1:K], lty=2, col=4)
lines(age, marg2.975[1:K], lty=2, col=4)
lines(age, marg2.[K+1:K], lwd=2, col=2)
lines(age, marg2.025[K+1:K], lty=2, col=2)
lines(age, marg2.975[K+1:K], lty=2, col=2)
dev.off()
##dev.copy2pdf(file='sort-growth.pdf')

                     
