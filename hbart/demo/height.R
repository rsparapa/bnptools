## marginal effects example
options(mc.cores=8)
library(hbart)

data(bmx)
bmx$RIDRETH2=factor(bmx$RIDRETH2)
str(bmx)
x.train=bmx[ , -(1:2)]
print(cor(bmx$BMXHT, x.train[ , c(2, 4)])^2)

ntree=300
ntreeh=60
numcut=200

file.='fit2-growth.rds'
if(file.exists(file.)) {
    XPtr=FALSE
    fit2=readRDS(file.)
} else {
    XPtr=TRUE
    set.seed(21)
    X.train=hbartModelMatrix(x.train)
    fit2 = hbart(X.train, bmx$BMXHT, ntree=ntree, ntreeh=ntreeh, numcut=numcut)
    saveRDS(fit2, file.)
}
print(sort(fit2$mu.varprob, TRUE))
print(cor(bmx$BMXHT, fit2$pred$f.train.mean)^2)

(N=length(bmx$BMXHT))
K=64
(age=seq(2, 18, length.out=K+1)[1:K])

x.train.=bmx[ , 3:4]
print(cor(bmx$BMXWT, x.train.[ , 2])^2)

file.='fit0-growth.rds'
if(file.exists(file.)) {
    XPtr=FALSE
    fit0=readRDS(file.)
} else {
    XPtr=TRUE
    set.seed(20)
    X.train. = hbartModelMatrix(x.train.)
    fit0 = hbart(X.train., bmx$BMXWT, ntree=ntree, ntreeh=ntreeh, numcut=numcut)
    saveRDS(fit0, file.)
}
print(cor(bmx$BMXWT, fit0$pred$f.train.mean)^2)

x.test. = x.train.
for(k in 2:K) x.test.=rbind(x.test., x.train.)
x.test.$RIDAGEEX=rep(age, each=N)
x.test.=rbind(x.test., x.test.)
x.test.$RIAGENDR=rep(1:2, each=N*K)
str(x.test.)

file.='pred0-growth.rds'
if(file.exists(file.)) {
    pred0=readRDS(file.)
} else {
    X.test. = hbartModelMatrix(x.test.)
    pred0 = predict(fit0, X.test., XPtr=XPtr)
    saveRDS(pred0, file.)
}

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

file.='pred2-growth.rds'
if(file.exists(file.)) {
    pred2=readRDS(file.)
} else {
    X.test0 =hbartModelMatrix(x.test0)
    pred2 = predict(fit2, X.test0, XPtr=XPtr)
    saveRDS(pred2, file.)
}

pred2.= pred2$f.test
prsd2.= pred2$s.test
marg2. = 0
marg2.025 = 0
marg2.10 = 0
marg2.25 = 0
marg2.75 = 0
marg2.90 = 0
marg2.975 = 0
for(k in 1:(2*K)) {
    marg2.[k] = mean(apply(pred2.[ , (k-1)*N+1:N], 1, mean))
    fpd = mean(apply(prsd2.[ , (k-1)*N+1:N], 1, mean))
    marg2.025[k] = marg2.[k]+qnorm(0.025)*fpd
    marg2.10[k] = marg2.[k]+qnorm(0.10)*fpd
    marg2.25[k] = marg2.[k]+qnorm(0.25)*fpd
    marg2.75[k] = marg2.[k]+qnorm(0.75)*fpd
    marg2.90[k] = marg2.[k]+qnorm(0.90)*fpd
    marg2.975[k] = marg2.[k]+qnorm(0.975)*fpd
}

col.=c(4, 2) ## males=blue, females=red

M = (bmx$RIAGENDR==1)
F = (bmx$RIAGENDR==2)

##pdf('M-growth.pdf')
plot(bmx$RIDAGEEX[M], bmx$BMXHT[M],
     pch='.', col=col.[bmx$RIAGENDR[M]],
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age, marg2.[1:K], lwd=2, col=4)
lines(age, marg2.025[1:K], lty=2, col=4)
lines(age, marg2.10[1:K], lty=2, col=4)
lines(age, marg2.25[1:K], lty=2, col=4)
lines(age, marg2.75[1:K], lty=2, col=4)
lines(age, marg2.90[1:K], lty=2, col=4)
lines(age, marg2.975[1:K], lty=2, col=4)
text(18, marg2.025[K], '2.5')
text(18, marg2.10[K], '10')
text(18, marg2.25[K], '25')
text(18, marg2.[K], '50')
text(18, marg2.75[K], '75')
text(18, marg2.90[K], '90')
text(18, marg2.975[K], '97.5')
##dev.off()

##pdf('F-growth.pdf')
plot(bmx$RIDAGEEX[F], bmx$BMXHT[F],
     pch='.', col=col.[bmx$RIAGENDR[F]],
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age, marg2.[K+1:K], lwd=2, col=2)
lines(age, marg2.025[K+1:K], lty=2, col=2)
lines(age, marg2.10[K+1:K], lty=2, col=2)
lines(age, marg2.25[K+1:K], lty=2, col=2)
lines(age, marg2.75[K+1:K], lty=2, col=2)
lines(age, marg2.90[K+1:K], lty=2, col=2)
lines(age, marg2.975[K+1:K], lty=2, col=2)
text(18, marg2.025[2*K], '2.5')
text(18, marg2.10[2*K], '10')
text(18, marg2.25[2*K], '25')
text(18, marg2.[2*K], '50')
text(18, marg2.75[2*K], '75')
text(18, marg2.90[2*K], '90')
text(18, marg2.975[2*K], '97.5')
##dev.off()

data(CDCheight)
M.=(CDCheight$sex==1 & CDCheight$age<18)
F.=(CDCheight$sex==2 & CDCheight$age<18)

pdf('M-CDCheight.pdf')
plot(bmx$RIDAGEEX[M], bmx$BMXHT[M], type='n',
     ylab='Height (cm)', sub='Boys',
     xlab='Age (yr)')
lines(age, marg2.[1:K], lwd=2, col=4)
lines(age, marg2.025[1:K], lty=2, col=4)
lines(age, marg2.10[1:K], lty=2, col=4)
lines(age, marg2.25[1:K], lty=2, col=4)
lines(age, marg2.75[1:K], lty=2, col=4)
lines(age, marg2.90[1:K], lty=2, col=4)
lines(age, marg2.975[1:K], lty=2, col=4)
lines(CDCheight$age[M.], CDCheight$height.50[M.], lwd=2)
lines(CDCheight$age[M.], CDCheight$height.025[M.], lty=2)
lines(CDCheight$age[M.], CDCheight$height.10[M.], lty=2)
lines(CDCheight$age[M.], CDCheight$height.25[M.], lty=2)
lines(CDCheight$age[M.], CDCheight$height.75[M.], lty=2)
lines(CDCheight$age[M.], CDCheight$height.90[M.], lty=2)
lines(CDCheight$age[M.], CDCheight$height.975[M.], lty=2)
text(18, marg2.025[K], '2.5')
text(18, marg2.10[K], '10')
text(18, marg2.25[K], '25')
text(18, marg2.[K], '50')
text(18, marg2.75[K], '75')
text(18, marg2.90[K], '90')
text(18, marg2.975[K], '97.5')
dev.off()

pdf('F-CDCheight.pdf')
plot(bmx$RIDAGEEX[F], bmx$BMXHT[F], type='n',
     ylab='Height (cm)', sub='Girls',
     xlab='Age (yr)')
lines(age, marg2.[K+1:K], lwd=2, col=2)
lines(age, marg2.025[K+1:K], lty=2, col=2)
lines(age, marg2.10[K+1:K], lty=2, col=2)
lines(age, marg2.25[K+1:K], lty=2, col=2)
lines(age, marg2.75[K+1:K], lty=2, col=2)
lines(age, marg2.90[K+1:K], lty=2, col=2)
lines(age, marg2.975[K+1:K], lty=2, col=2)
lines(CDCheight$age[F.], CDCheight$height.50[F.], lwd=2)
lines(CDCheight$age[F.], CDCheight$height.025[F.], lty=2)
lines(CDCheight$age[F.], CDCheight$height.10[F.], lty=2)
lines(CDCheight$age[F.], CDCheight$height.25[F.], lty=2)
lines(CDCheight$age[F.], CDCheight$height.75[F.], lty=2)
lines(CDCheight$age[F.], CDCheight$height.90[F.], lty=2)
lines(CDCheight$age[F.], CDCheight$height.975[F.], lty=2)
text(18, marg2.025[2*K], '2.5')
text(18, marg2.10[2*K], '10')
text(18, marg2.25[2*K], '25')
text(18, marg2.[2*K], '50')
text(18, marg2.75[2*K], '75')
text(18, marg2.90[2*K], '90')
text(18, marg2.975[2*K], '97.5')
dev.off()

