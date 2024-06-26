## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)
str(bmx)

strata <- stratrs(bmx$BMXHT, 5, seed = 10)
tapply(bmx$BMXHT, strata, mean) ## height

tapply(bmx$BMXWT, strata, mean) ## weight
tapply(bmx$RIDAGEEX, strata, mean) ## age
tapply(bmx$RIAGENDR, strata, mean) ## 1:M, 2:F

train <- (strata<5)
x.train=bmx[train, -(1:2)]
str(x.train)
x.test =bmx[!train, -(1:2)]
dim(x.test)

print(cor(bmx$BMXHT[train], x.train[ , 2:3])^2)

file.='growth1-fit1.rds'
if(!file.exists(file.)) {
    ## one long chain: single thread computing
    ## set.seed(21)
    ## fit1 = gbart(x.train, bmx$BMXHT[train], x.test=x.test) 
    ## multiple shorter chains: parallel multi-threading
    fit1 = mc.gbart(x.train, bmx$BMXHT[train], x.test=x.test, seed=21)
    saveRDS(fit1, file.)
} else { fit1=readRDS(file.) }
print(c(train.Rsqr = cor(bmx$BMXHT[train], fit1$yhat.train.mean)^2))
print(c(test.Rsqr = cor(bmx$BMXHT[!train], fit1$yhat.test.mean)^2))

col.=c(4, 2) ## 1:M=blue, 2:F=red
plot(fit1$yhat.test.mean, bmx$BMXHT[!train], asp=1,
     pch=19, col=col.[bmx$RIAGENDR[!train]],
     cex = log2(bmx$RIDAGEEX[!train])/2,
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
points(fit1$yhat.test.mean, bmx$BMXHT[!train], 
     cex = log2(bmx$RIDAGEEX[!train])/2, pch=21)
abline(b=1, a=0, col=8)
age. <- seq(2, 18, 4)
legend('bottomright', paste0(age.),
       pch = 1, cex = log2(age.)/2, horiz = TRUE)
legend('bottomleft', c(' ', 'M', 'F'), col = c(0, col.), pch = 19)
##dev.copy2pdf(file='growth1-fit1.pdf')

age. <- 2:17
x.test. <- cbind(0, age.) ## column 1: sex, column 2: age
K <- nrow(x.test.)
x.test. <- rbind(x.test., x.test.)
x.test.[ , 1] <- rep(1:2, each = K)
x.test.

print(system.time(pred1 <- FPD(fit1, x.test., S = 1:2))) ## assuming independence

F.test <- (!train & bmx$RIAGENDR == 2)
plot(bmx$RIDAGEEX[F.test], bmx$BMXHT[F.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred1$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred1$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred1$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-indF.pdf')

M.test <- (!train & bmx$RIAGENDR == 1)
plot(bmx$RIDAGEEX[M.test], bmx$BMXHT[M.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age.,  pred1$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred1$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred1$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-indM.pdf')

file.='growth1-fit2.rds'
if(!file.exists(file.)) {
    ## set.seed(22)
    ## fit2 = gbart(x.train[ , 1:2], bmx$BMXWT[train], x.test=x.test.) 
    fit2 = mc.gbart(x.train[ , 1:2], bmx$BMXWT[train], x.test=x.test., seed=22)
    saveRDS(fit2, file.)
} else { fit2=readRDS(file.) }
print(c(train.Rsqr = cor(bmx$BMXWT[train], fit2$yhat.train.mean)^2))

(x.test. <- cbind(x.test., fit2$yhat.test.mean))
print(system.time(pred2 <- FPD(fit1, x.test., S = 1:3))) ## conditional means

plot(bmx$RIDAGEEX[F.test], bmx$BMXHT[F.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred2$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmF.pdf')

plot(bmx$RIDAGEEX[M.test], bmx$BMXHT[M.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred2$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmM.pdf')

(x.test. <- x.test.[ , -3])
print(system.time(pred3 <- FPD(fit1, x.test., S = 1:2, subset. = 1:2))) ## conditional dependence

plot(bmx$RIDAGEEX[F.test], bmx$BMXHT[F.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred3$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred3$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred3$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cdF.pdf')

plot(bmx$RIDAGEEX[M.test], bmx$BMXHT[M.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred3$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred3$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred3$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cdM.pdf')

