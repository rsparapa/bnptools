
options(mc.cores=8)
library(nftbart)

data(bmx)
str(bmx)

x.train=bmx[ , 3:4]
age=seq(2, 17.75, 0.25)
K=length(age)
x.test=rbind(cbind(1, age), cbind(2, age))
dimnames(x.test)[[2]]=names(x.train)
set.seed(21)
post=nft(x.train, times=bmx$BMXHT, x.test=x.test)

print(summary(bmx$BMXHT))
print(cor(bmx$BMXHT, exp(post$f.train.mean))^2)

Q=c(2.5, 10, 25, 50, 75, 90, 97.5)
L=length(Q)
H=ncol(post$dpwt.)
q=matrix(nrow=post$ndpost, ncol=L)
for(j in 1:L) {
    for(k in 1:H)
        if(k==1)
            q[ , j]=post$dpwt.[ , k]*
                qnorm(Q[j]/100, post$dpmu.[ , k], post$dpsd.[ , k])
    else q[ , j]=q[ , j]+
             post$dpwt.[ , k]*
             qnorm(Q[j]/100, post$dpmu.[ , k], post$dpsd.[ , k])
    }
    
q=apply(q, 2, mean)    
##q=c(-1.96, -1.28, -0.67, 0, 0.67, 1.28, 1.96)
lty=c(2, 2, 2, 1, 2, 2, 2)
col=c(4, 2)
M=(bmx$RIAGENDR==1)
F=(!M)

CDC = read.csv('height.csv')                     
str(CDC)
CDC$age=CDC$month/12
## sex coded like EPIC: the opposite of NHANES
M.=(CDC$sex_c==2 & CDC$age<18)
F.=(CDC$sex_c==1 & CDC$age<18)

pdf(file='bmx-M.pdf')
plot(bmx$RIDAGEEX[M], bmx$BMXHT[M], pch='.', col=4,
     xlab='Age', ylab='Height (cm)')
for(j in 1:L) {
    lines(age, sort(exp(post$f.test.mean[1:K]+q[j]*post$s.test.mean[1:K])),
          col=4, lty=lty[j], lwd=3-lty[j])
    text(age[K]+0.25, exp(post$f.test.mean[K]+q[j]*post$s.test.mean[K]), Q[j])
}
dev.off()

pdf(file='bmx-F.pdf')
plot(bmx$RIDAGEEX[F], bmx$BMXHT[F], pch='.', col=2,
     xlab='Age', ylab='Height (cm)')
for(j in 1:L) {
    lines(age, sort(exp(post$f.test.mean[K+1:K]+q[j]*post$s.test.mean[K+1:K])),
          col=2, lty=lty[j], lwd=3-lty[j])
    text(age[K]+0.25, exp(post$f.test.mean[2*K]+q[j]*post$s.test.mean[2*K]), Q[j])
}
dev.off()
