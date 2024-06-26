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
print(system.time(pred2 <- FPD(fit1, x.test., S = 1:3))) ## conditional means

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

file. <- "growth3-pred.t.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.t <- SHAP(fit1, x.test.[1:L, 2], S = 2)))
    saveRDS(pred.t, file.) 
} else {
    pred.t <- readRDS(file.) 
}

pred.t$yhat.test.mean <- pred.t$yhat.test.mean+fit1$offset

file. <- "growth3-pred.w.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.w <- SHAP(fit1, x.test.[ , 3], S = 3)))
    saveRDS(pred.w, "growth3-pred.w.RDS")
} else {
    pred.w <- readRDS(file.) 
}

file. <- "growth3-pred.u.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.u <- SHAP(fit1, 1:2, S = 1)))
    saveRDS(pred.u, "growth3-pred.u.RDS")
} else {
    pred.u <- readRDS(file.) 
}

plot(bmx$RIDAGEEX, bmx$BMXHT, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[L+1:L], col = 2, lwd = 2)
lines(age., pred2$yhat.test.mean[1:L], col = 4, lwd = 2)
points(age., pred.t$yhat.test.mean)
points(age., pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[L+1:L], col = 2, pch = 25)
points(age., pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[1:L], col = 4, pch = 24)
legend('topleft', legend = c('FPD Conditional Means', 'SHAP: age-only', 
'SHAP: age, F, weight', 'SHAP: age, M, weight'),
       pch = c(32, 1, 25, 24), col = c(0, 1, 2, 4))
dev.copy2pdf(file='growth3-SHAP.pdf')

