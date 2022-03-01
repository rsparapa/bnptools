## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

L=1000 ## 4000 burn-in discard
T=10   ## thin 
M=1000 ## 2000 keep

file.='height3-fit2.rds'
if(file.exists(file.)) { fit2=readRDS(file.)
} else {
    fit2 = mc.gbart(x.train, bmx$BMXHT, seed=21, sparse=TRUE,
                    ndpost=M, nskip=L, keepevery=T)
    saveRDS(fit2, file.)
}
print(sort(fit2$varprob.mean, TRUE))

##source('Rhat.R')
## check=maxRhat(c(fit2$sigma), fit2$chains)
## print(check$maxRhat)
## plot(check$rho[1:50], type='h', ylim=c(-1, 1))
## plot(acf(c(fit2$sigma)), ylim=c(-1, 1))

pdf('height3-zep.pdf')
col.=c(4, 2) ## males=blue, females=red
print(cor(bmx$BMXHT, fit2$yhat.train.mean)^2)
plot(fit2$yhat.train.mean, bmx$BMXHT, asp=1,
     pch='.', col=col.[bmx$RIAGENDR],
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
##legend('topleft', col=c(4, 2), legend=c('M', 'F'))
text(80, 160, expression(italic(R)^2), pos=4)
text(85, 160, '=96.7%', pos=4)
abline(b=1, a=0, col=8)
dev.off()

pdf('height3-sigma.pdf')
plot(fit2$sigma[ , 1], type='l', ##log='y',
     ylim=c(0, 10),
     ylab=expression(italic(sigma)),
     xlab=paste0('MCMC sequence: ', fit2$chains, ' chain(s)'))
for(i in 2:8) lines(fit2$sigma[ , i], col=i)
abline(v=c(0, L), h=c(0, fit2$sigma.mean), lty=2)
dev.off()

(N=length(bmx$BMXHT))
K=64
(age=seq(2, 18, length.out=K+1)[1:K])

x.train.=bmx[ , 3:4]
print(cor(bmx$BMXWT, x.train.[ , 2])^2)

file.='height3-fit0.rds'
if(file.exists(file.)) { fit0=readRDS(file.)
} else {
    fit0 = mc.gbart(x.train., bmx$BMXWT, seed=20, sparse=TRUE,
                    ndpost=M, nskip=L, keepevery=T)
    saveRDS(fit0, file.)
}
print(cor(bmx$BMXWT, fit0$yhat.train.mean)^2)

x.test. = x.train.
for(k in 2:K) x.test.=rbind(x.test., x.train.)
x.test.$RIDAGEEX=rep(age, each=N)
x.test.=rbind(x.test., x.test.)
x.test.$RIAGENDR=rep(1:2, each=N*K)
str(x.test.)

file.='height3-pred0.rds'
if(file.exists(file.)) { pred0=readRDS(file.)
} else {
    pred0 = predict(fit0, x.test.)
    saveRDS(pred0, file.)
}

(M = nrow(pred0))
(L = ncol(pred0)/(N*K))
marg0 = 0
for(l in 1:L) { 
    h=(l-1)*K
    for(k in 1:K) {
        j=h+k
        marg0[j] = mean(apply(pred0[ , (j-1)*N+1:N], 1, mean))
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

file.='height3-pred2.rds'
if(file.exists(file.)) { pred2.=readRDS(file.)
} else {
    pred2. = predict(fit2, x.test0)
    saveRDS(pred2., file.)
}

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
plot(bmx$RIDAGEEX,  bmx$BMXHT,
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

                     
