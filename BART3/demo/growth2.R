## marginal effects example
options(mc.cores=8)
library(BART3)

data(bmx)

strata <- stratrs(bmx$BMXHT, 5, seed = 10)
train <- (strata<5)
x.train=bmx[train, -(1:2)]
x.test =bmx[!train, -(1:2)]

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

age. <- 2:17
x.test. <- cbind(0, age.) ## column 1: sex, column 2: age
L <- nrow(x.test.)
x.test. <- rbind(x.test., x.test.)
x.test.[ , 1] <- rep(1:2, each = L)

file.='growth1-fit2.rds'
if(!file.exists(file.)) {
    ## set.seed(22)
    ## fit2 = gbart(x.train[ , 1:2], bmx$BMXWT[train], x.test=x.test.) 
    fit2 = mc.gbart(x.train[ , 1:2], bmx$BMXWT[train], x.test=x.test., seed=22)
    saveRDS(fit2, file.)
} else { fit2=readRDS(file.) }
print(c(train.Rsqr = cor(bmx$BMXWT[train], fit2$yhat.train.mean)^2))

(x.test. <- cbind(x.test., fit2$yhat.test.mean))
print(system.time(pred2 <- FPD(fit1, x.test., S = 1:3))) ## Imputation Marginal

F.test <- (!train & bmx$RIAGENDR == 2)
plot(bmx$RIDAGEEX[F.test], bmx$BMXHT[F.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[L+1:L], col=2, lwd=2)
lines(age., pred2$yhat.test.lower[L+1:L], col=2, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[L+1:L], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmF.pdf')

M.test <- (!train & bmx$RIAGENDR == 1)
plot(bmx$RIDAGEEX[M.test], bmx$BMXHT[M.test],
     pch=1, col=8,
     cex = log2(bmx$RIDAGEEX[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[1:L], col=4, lwd=2)
lines(age., pred2$yhat.test.lower[1:L], col=4, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[1:L], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmM.pdf')

print(system.time(pred30 <- FPDK(fit1, x.test., S = 1:3, seed = 22))) ## K=30
print(system.time(pred5 <- FPDK(fit1, x.test., S = 1:3, seed = 22, mult.impute = 5))) ## K=5

plot(bmx$RIDAGEEX, bmx$BMXHT, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)')
lines(c(1.65, 1.95), c(188, 188), lwd = 2)
lines(age., pred2$yhat.test.mean[L+1:L], col = 2, lwd = 2)
lines(age., pred2$yhat.test.mean[1:L], col = 4, lwd = 2)
points(age., pred30$yhat.test.mean[L+1:L], col = 2, pch = 24)
points(age., pred30$yhat.test.mean[1:L], col = 4, pch = 24)
points(age., pred5$yhat.test.mean[L+1:L], col = 2, pch = 25)
points(age., pred5$yhat.test.mean[1:L], col = 4, pch = 25)
legend('topleft', legend = c('Imputation Marginal', 'FPD: N=2748', 'FPDK: K=30', 'FPDK: K=5'),
       pch = c(32, 32, 24, 25))
##dev.copy2pdf(file='growth2-FPDK.pdf')

plot(bmx$RIDAGEEX, bmx$BMXHT, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)',
     xlim = c(12, 17), ylim = c(150, 180))
lines(c(11.875, 12), c(178.75, 178.75), lwd = 2)
lines(age., pred2$yhat.test.lower[L+1:L], col = 2, lwd = 2)
lines(age., pred2$yhat.test.upper[L+1:L], col = 2, lwd = 2)
lines(age., pred2$yhat.test.lower[1:L], col = 4, lwd = 2)
lines(age., pred2$yhat.test.upper[1:L], col = 4, lwd = 2)
points(age., pred30$yhat.test.upper[L+1:L], col = 2, pch = 24)
points(age., pred30$yhat.test.lower[L+1:L], col = 2, pch = 24)
points(age., pred30$yhat.test.upper[1:L], col = 4, pch = 24)
points(age., pred30$yhat.test.lower[1:L], col = 4, pch = 24)
points(age., pred5$yhat.test.upper[L+1:L], col = 2, pch = 25)
points(age., pred5$yhat.test.lower[L+1:L], col = 2, pch = 25)
points(age., pred5$yhat.test.upper[1:L], col = 4, pch = 25)
points(age., pred5$yhat.test.lower[1:L], col = 4, pch = 25)
legend('topleft', legend = c('Imputation Marginal', 
                             'FPD: N=2748', 'FPDK: K=30', 'FPDK: K=5'),
       pch = c(32, 32, 24, 25))
##dev.copy2pdf(file='growth2-FPDK-EV.pdf')
