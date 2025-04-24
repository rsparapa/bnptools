## marginal effects example
library(lattice)
options(mc.cores=8)
library(BART3)

data(bmx)
str(bmx)

data(CDCweight)
##str(CDCweight)

SDwt <- tapply(CDCweight$SD, CDCweight$age, mean)
age <- as.double(names(SDwt))
SDwt <- SDwt/SDwt[age == 8]
plot(age, SDwt, type = 'l', lwd = 2, xlim = c(2, 18), log = 'y', lty = 2,
     ylim = c(1/3, 3), ylab = 'relative uncertainty')
abline(h = 1, col = 8)

SDwt <- data.frame(SDwt = SDwt, months = round(age*12))
str(SDwt)

data(bmx)
bmx$months <- round(bmx$age*12)

bmxwt <- merge(bmx, SDwt, by = 'months')

data(CDCheight)

SDht <- tapply(CDCheight$SD, CDCheight$age, mean)
age <- as.double(names(SDht))
SDht <- SDht/SDht[age == 8]
lines(age, SDht, lty = 1, lwd = 2)
legend('topleft', c('height', 'weight'), lwd = 2, lty = 1:2)

SDht <- data.frame(SDht = SDht, months = round(age*12))
str(SDht)

bmxht <- merge(bmxwt, SDht, by = 'months')

(N <- nrow(bmxht))

set.seed(10)
strata <- rep(0, N)
M <- (bmxht$sex == 1)
F <- (bmxht$sex == 2)
strata.M <- stratrs(bmxht$height[M], 5)
strata.F <- stratrs(bmxht$height[F], 5)
strata[M] <- strata.M
strata[F] <- strata.F
table(strata, M)

train <- (strata<5)
pick <- c('age', 'sex', 'weight', 'race') ## t, u, w, v
x.train <- bmxht[train, pick]
x.test  <- bmxht[!train, pick]

file.='child1-fit1.rds'
if(!file.exists(file.)) {
    ## one long chain: single thread computing
    ## set.seed(21)
    ## fit1 = gbart(x.train, bmx$BMXHT[train], x.test=x.test) 
    ## multiple shorter chains: parallel multi-threading
    fit1 = mc.gbart(x.train, bmxht$height[train], x.test=x.test, 
                    seed=21, w = bmxht$SDht[train])
    saveRDS(fit1, file.)
} else { fit1=readRDS(file.) }
print(c(train.Rsqr = cor(bmxht$height[train], fit1$yhat.train.mean)^2)) ## 96.6
print(c(test.Rsqr = cor(bmxht$height[!train], fit1$yhat.test.mean)^2))  ## 96.3

age. <- 2:17
K <- length(age.)
x.test. <- cbind(age., 0) 
x.test. <- rbind(x.test., x.test.)
x.test.[ , 2] <- rep(1:2, each = K)

##print(system.time(pred1 <- FPD(fit1, x.test., S = 1:2))) ## assuming independence

file.='child1-fit2.rds'
if(!file.exists(file.)) {
    ## set.seed(22)
    ## fit2 = gbart(x.train[ , 1:2], bmx$BMXWT[train], x.test=x.test.) 
    fit2 = mc.gbart(x.train[ , 1:2], bmxht$weight[train], 
                    x.test=x.test.[ , 1:2], seed=22, w = bmxht$SDwt[train])
    saveRDS(fit2, file.)
} else { fit2=readRDS(file.) }
print(c(train.Rsqr = cor(bmxht$weight[train], fit2$yhat.train.mean)^2)) ## 71.3

(x.test. <- cbind(x.test.[ , 1:2], fit2$yhat.test.mean))

print(system.time(pred2 <- FPD(fit1, x.test., S = c('age', 'sex', 'weight'))))  ## synthetic approx

F.test <- (!train & bmxht$sex == 2)
plot(bmxht$age[F.test], bmxht$height[F.test],
     pch=1, col=8,
     cex = log2(bmx$age[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred2$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cmF.pdf')

M.test <- (!train & bmxht$sex == 1)
plot(bmxht$age[M.test], bmxht$height[M.test],
     pch=1, col=8,
     cex = log2(bmxht$age[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred2$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cmM.pdf')

file. <- "child2-pred.tu.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.tu <- SHNN2(fit1, x.test.[ , 1:2], S = c('age', 'sex'))))
    saveRDS(pred.tu, file.) 
} else {
    pred.tu <- readRDS(file.) 
}
##     user   system  elapsed 
## 4862.598    2.152 4864.989 

pred.tu. <- list()
pred.tu.$yhat.test.mean <- apply(pred.tu, 2, mean)
pred.tu.$yhat.test.lower <- apply(pred.tu, 2, quantile, probs = 0.025)
pred.tu.$yhat.test.upper <- apply(pred.tu, 2, quantile, probs = 0.975)

file. <- "child2-pred.t.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.t <- SHNN(fit1, age., S = 'age')))
    saveRDS(pred.t, file.) 
} else {
    pred.t <- readRDS(file.) 
}

pred.t$yhat.test <- pred.t$yhat.test+fit1$offset
pred.t$yhat.test.mean <- pred.t$yhat.test.mean+fit1$offset
pred.t$yhat.test.lower <- pred.t$yhat.test.lower+fit1$offset
pred.t$yhat.test.upper <- pred.t$yhat.test.upper+fit1$offset

file. <- "child2-pred.w.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.w <- SHNN(fit1, x.test.[ , 3], S = 'weight')))
    saveRDS(pred.w, file.)
} else {
    pred.w <- readRDS(file.) 
}

file. <- "child2-pred.u.RDS" 
if(!file.exists(file.)) {
    print(system.time(pred.u <- SHNN(fit1, 1:2, S = 'sex')))
    saveRDS(pred.u, file.)
} else {
    pred.u <- readRDS(file.) 
}

##pdf(file = 'growth6.pdf')
plot(bmx$age, bmx$height, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col = 2, lwd = 2)
lines(age., pred2$yhat.test.mean[1:K], col = 4, lwd = 2)
points(age., pred.t$yhat.test.mean)
points(age., pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[K+1:K], col = 2, pch = 25)
points(age., pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[1:K], col = 4, pch = 24)
lines(c(1.65, 1.95), c(188, 188), lwd = 2)
legend('topleft', legend = c('Imputation Marginal', 'FPD', 'SHAP: age-only', 
'SHAP: age, F, weight', 'SHAP: age, M, weight'),
       pch = c(32, 32, 1, 25, 24), col = c(0, 0, 1, 2, 4))
##dev.copy2pdf(file='growth3-SHAP.pdf')

plot(bmx$age, bmx$height, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col = 2, lwd = 2)
lines(age., pred2$yhat.test.mean[1:K], col = 4, lwd = 2)
points(age., pred.t$yhat.test.mean)
points(age., 2*pred.tu[K+1:K]+pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[K+1:K], col = 2, pch = 25)
points(age., 2*pred.tu[1:K]+pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[1:K], col = 4, pch = 24)
lines(c(1.65, 1.95), c(188, 188), lwd = 2)
legend('topleft', legend = c('Imputation Marginal', 'FPD', 'SHAP: age-only', 
'SHAP: age, F, weight', 'SHAP: age, M, weight'),
       pch = c(32, 32, 1, 25, 24), col = c(0, 0, 1, 2, 4))
##dev.copy2pdf(file='growth3-SHAP.pdf')
##dev.off()

plot(bmx$age, bmx$height, type = 'n',
     ylab='Height (cm)', xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col = 2, lwd = 2)
lines(age., pred2$yhat.test.mean[1:K], col = 4, lwd = 2)
points(age., pred.t$yhat.test.mean)
points(age., 2*(pred.tu[K+1:K])+pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[K+1:K], col = 2, pch = 25)
points(age., 2*(pred.tu[1:K]  )+pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[1:K], col = 4, pch = 24)
lines(c(1.65, 1.95), c(188, 188), lwd = 2)
legend('topleft', legend = c('Imputation Marginal', 'FPD', 'SHAP: age-only', 
'SHAP: age, F, weight', 'SHAP: age, M, weight'),
       pch = c(32, 32, 1, 25, 24), col = c(0, 0, 1, 2, 4))
##dev.copy2pdf(file='growth3-SHAP.pdf')

df <- data.frame(sex = bmx$sex[train], age = bmx$age[train], height = bmx$height[train])
N <- nrow(df)
df <- rbind(df, df, df)
df$marginal <- rep(c('A', 'DI', 'NN'), each = N)
df$marginal <- factor(df$marginal)

pdf(file = 'SHAP.pdf')
xyplot(height~age|sex+marginal, data = df, as.table = TRUE, strip = FALSE, col = 8, type = 'n', ylab = 'height (cm)',
xlim = c(1.5, 17.5),
panel=function(...) { 
           i=panel.number()
           M <- (df$sex == 1 & df$marginal == 'A')
           F <- (df$sex == 2 & df$marginal == 'A')
           black. <- adjustcolor("black",alpha.f=0.2)
           grey. <- adjustcolor("gray",alpha.f=0.2)
           if((i%%2) == 1) {
               lpoints(df$age[M], df$height[M], col = grey., pch = 21)
           } else if((i%%2) == 0) {
               lpoints(df$age[F], df$height[F], col = grey., pch = 21)
           }
           if(i == 1) { 
               ##llines(age., pred2$yhat.test.mean[1:K], col = 1)
               llines(age., pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean, col = 4, lwd = 2)
               llines(age., apply(pred.u$yhat.test[ , 1]+pred.t$yhat.test, 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(pred.u$yhat.test[ , 1]+pred.t$yhat.test, 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
               ltext(2, 180, 'SV: t,u', pos = 4)
               ltext(15, 85, 'M', col = 4, pos = 4)
           } else if(i == 2) { 
               ##llines(age., pred2$yhat.test.mean[K+1:K], col = 1)
               llines(age., pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean, col = 2, lwd = 2)
               ##llines(age., pred.t$yhat.test.lower, col = 2, lty = 3, lwd = 2)
               ##llines(age., pred.t$yhat.test.upper, col = 2, lty = 3, lwd = 2)
               llines(age., apply(pred.u$yhat.test[ , 2]+pred.t$yhat.test, 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(pred.u$yhat.test[ , 2]+pred.t$yhat.test, 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
               ltext(15, 85, 'F', col = 2, pos = 4)
           } else if(i == 3) { 
               ##llines(age., pred2$yhat.test.mean[1:K], col = 1)
               ltext(2, 180, 'SV: t,u,t:u', pos = 4)
               llines(age., 2*pred.tu.$yhat.test.mean[1:K]+pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean, col = 4, lwd = 2)
               llines(age., apply(2*pred.tu[ , 1:K]+pred.u$yhat.test[ , 1]+pred.t$yhat.test, 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*pred.tu[ , 1:K]+pred.u$yhat.test[ , 1]+pred.t$yhat.test, 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 4) { 
               ##llines(age., pred2$yhat.test.mean[K+1:K], col = 1)
               llines(age., 2*pred.tu.$yhat.test.mean[K+1:K]+pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean, col = 2, lwd = 2)
               llines(age., apply(2*pred.tu[ , K+1:K]+pred.u$yhat.test[ , 2]+pred.t$yhat.test, 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*pred.tu[ , K+1:K]+pred.u$yhat.test[ , 2]+pred.t$yhat.test, 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           } else if(i == 5) { 
               ##llines(age., pred2$yhat.test.mean[1:K], col = 1)
               ltext(2, 180, 'SA: t,u,t:u,w(t,u)', pos = 4)
               llines(age., 2*pred.tu.$yhat.test.mean[1:K]+pred.u$yhat.test.mean[1]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., apply(2*pred.tu[ , 1:K]+pred.u$yhat.test[ , 1]+pred.t$yhat.test+pred.w$yhat.test[, 1:K], 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*pred.tu[ , 1:K]+pred.u$yhat.test[ , 1]+pred.t$yhat.test+pred.w$yhat.test[ , 1:K], 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 6) { 
               ##llines(age., pred2$yhat.test.mean[K+1:K], col = 1)
               llines(age., 2*pred.tu.$yhat.test.mean[K+1:K]+pred.u$yhat.test.mean[2]+pred.t$yhat.test.mean+pred.w$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., apply(2*pred.tu[ , K+1:K]+pred.u$yhat.test[ , 2]+pred.t$yhat.test+pred.w$yhat.test[, K+1:K], 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*pred.tu[ , K+1:K]+pred.u$yhat.test[ , 2]+pred.t$yhat.test+pred.w$yhat.test[ , K+1:K], 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           }
})
dev.off()
     
