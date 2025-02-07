## marginal effects 
options(mc.cores = 8)
library(BART3)
library(lattice)

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

col.=c(4, 2) ## 1:M=blue, 2:F=red
plot(fit1$yhat.test.mean, bmxht$height[!train], asp=1,
     pch=19, col=col.[bmxht$sex[!train]],
     cex = log2(bmxht$age[!train])/2,
     xlab='Predicted Height (cm)',
     ylab='Observed Height (cm)')
points(fit1$yhat.test.mean, bmxht$height[!train], 
     cex = log2(bmxht$age[!train])/2, pch=21)
abline(b=1, a=0, col=8)
age. <- seq(2, 18, 4)
legend('bottomright', paste0(age.),
       pch = 1, cex = log2(age.)/2, horiz = TRUE)
legend('bottomleft', c(' ', 'M', 'F'), col = c(0, col.), pch = 19)
##dev.copy2pdf(file='child1-fit1.pdf')

age. <- 2:17
x.test. <- cbind(age., 0) ## age, sex
K <- nrow(x.test.)
x.test. <- rbind(x.test., x.test.)
x.test.[ , 2] <- rep(1:2, each = K)
x.test.

print(system.time(pred1 <- FPD(fit1, x.test., S = c('age', 'sex')))) ## assuming independence
##   user  system elapsed 
##230.347   1.132  51.217 

F.test <- (!train & bmxht$sex == 2)
plot(bmxht$age[F.test], bmxht$height[F.test],
     pch=1, col=8,
     cex = log2(bmxht$age[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred1$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred1$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred1$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='child1-indF.pdf')

M.test <- (!train & bmxht$sex == 1)
plot(bmxht$age[M.test], bmxht$height[M.test],
     pch=1, col=8,
     cex = log2(bmxht$age[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred1$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred1$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred1$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='child1-indM.pdf')

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
##   user  system elapsed 
##220.346   1.195  49.972
print(system.time(pred4 <- SRS(fit1, x.test.[ , 1:2], ##kern.var = FALSE,
                  S = c('age', 'sex')))) 
                  ##S = c('age', 'sex', 'weight')))) 
##   user  system elapsed 
##158.179   2.209  20.529

plot(bmxht$age[F.test], bmxht$height[F.test],
     pch=1, col=8,
     cex = log2(bmx$age[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred2$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cmF.pdf')

plot(bmxht$age[M.test], bmxht$height[M.test],
     pch=1, col=8,
     cex = log2(bmxht$age[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred2$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred2$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred2$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cmM.pdf')

(x.test. <- x.test.[ , -3])
print(system.time(pred3 <- NN(fit1, x.test., S = c('age', 'sex'), 
                              nearest = c('age', 'sex'))))
##   user  system elapsed 
## 30.123   0.965  22.040

plot(bmxht$age[F.test], bmxht$height[F.test],
     pch=1, col=8,
     cex = log2(bmxht$age[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred3$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred3$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred3$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cdF.pdf')

plot(bmxht$age[M.test], bmxht$height[M.test],
     pch=1, col=8,
     cex = log2(bmxht$age[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred3$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred3$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred3$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='child1-cdM.pdf')

## by repeating pred M times: fast
M <- 30 
file.='child1-pred0.rds'
if(!file.exists(file.)) {
    relwt <- data.frame(age = age., months = 12*age.)
    relwt <- merge(relwt, SDht, by = 'months')
    relwt <- rep(relwt$SDht, 2*M)
    set.seed(10)
    i <- sample.int(nrow(x.train), 2*K*M)
    x.test. <- cbind(rep(x.test.[ , 1], M), rep(x.test.[ , 2], M),
                     rnorm(2*K*M, rep(fit2$yhat.test.mean, M), 
                           relwt*fit2$sigma.mean),
                     fit1$x.train[i, 4:6])
    pred0 <- predict(fit1, x.test.) 
    saveRDS(pred0, file.)
} else { pred0=readRDS(file.) }

h <- seq(1, 2*K*M, 2*K)
pred <- list(yhat.test = pred0)
pred0 <- list()
pred0$yhat.test <- cbind(apply(pred$yhat.test[ , h], 1, mean))
for(i in 2:(2*K))
    pred0$yhat.test <- cbind(pred0$yhat.test,
                             apply(pred$yhat.test[ , h+i-1], 1, mean))

    pred0$yhat.test.mean  <- apply(pred0$yhat.test, 2, mean)
    pred0$yhat.test.lower <- apply(pred0$yhat.test, 2, quantile, 0.025)
    pred0$yhat.test.upper <- apply(pred0$yhat.test, 2, quantile, 0.975)

## by repeating FPD M times: slow
## M <- 30 
## file.='child1-pred.rds'
## if(!file.exists(file.)) {
##     relwt <- data.frame(age = age., months = 12*age.)
##     relwt <- merge(relwt, SDht, by = 'months')
##     relwt <- rep(relwt$SDht, 2*M)
##     set.seed(10)
##     x.test. <- cbind(rep(x.test.[ , 1], M), rep(x.test.[ , 2], M),
##                      rnorm(2*K*M, rep(fit2$yhat.test.mean, M), 
##                            relwt*fit2$sigma.mean))
##                            ##relwt*fit2$sigma.mean/sqrt(sum(train))))
##     pred <- FPD(fit1, x.test., S = c('age', 'sex', 'weight')) 
##     saveRDS(pred, file.)
## } else { pred=readRDS(file.) }

## h <- seq(1, 2*K*M, 2*K)
## pred0 <- list()
## pred0$yhat.test <- cbind(apply(pred$yhat.test[ , h], 1, mean))
## for(i in 2:(2*K))
##     pred0$yhat.test <- cbind(pred0$yhat.test,
##                              apply(pred$yhat.test[ , h+i-1], 1, mean))

##     pred0$yhat.test.mean  <- apply(pred0$yhat.test, 2, mean)
##     pred0$yhat.test.lower <- apply(pred0$yhat.test, 2, quantile, 0.025)
##     pred0$yhat.test.upper <- apply(pred0$yhat.test, 2, quantile, 0.975)

K <- length(age.)
plot(bmxht$age[F.test], bmxht$height[F.test],
     pch=1, col=8,
     cex = log2(bmxht$age[F.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred0$yhat.test.mean[K+1:K], col=2, lwd=2)
lines(age., pred0$yhat.test.lower[K+1:K], col=2, lwd=2, lty=2)
lines(age., pred0$yhat.test.upper[K+1:K], col=2, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmF.pdf')

plot(bmxht$age[M.test], bmxht$height[M.test],
     pch=1, col=8,
     cex = log2(bmxht$age[M.test])/2,
     ylab='Height (cm)',
     xlab='Age (yr)')
lines(age., pred0$yhat.test.mean[1:K], col=4, lwd=2)
lines(age., pred0$yhat.test.lower[1:K], col=4, lwd=2, lty=2)
lines(age., pred0$yhat.test.upper[1:K], col=4, lwd=2, lty=2)
##dev.copy2pdf(file='growth1-cmM.pdf')


df <- data.frame(sex = bmxht$sex[train], age = bmxht$age[train], height = bmxht$height[train])
N <- nrow(df)
df <- rbind(df, df, df, df, df)
df$marginal <- rep(c('A', 'B', 'C', 'D', 'RS'), each = N)
df$marginal <- factor(df$marginal)

pdf(file = 'FPD-final.pdf')
xyplot(height~age|sex+marginal, data = df, as.table = TRUE, strip = FALSE,
xlim = c(1.5, 17.5),
       col = 8, type = 'n', ylab = 'height (cm)', panel=function(...) { 
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
               llines(age., pred1$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., pred1$yhat.test.lower[1:K], col = 4, lwd = 2, lty = 3)
               llines(age., pred1$yhat.test.upper[1:K], col = 4, lwd = 2, lty = 3)
               ltext(3, 180, 'FPD: t,u', pos = 4)
               ltext(15, 85, 'M', col = 4, pos = 4)
           } else if(i == 2) { 
               llines(age., pred1$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., pred1$yhat.test.lower[K+1:K], col = 2, lwd = 2, lty = 3)
               llines(age., pred1$yhat.test.upper[K+1:K], col = 2, lwd = 2, lty = 3)
               ltext(15, 85, 'F', col = 2, pos = 4)
           } else if(i == 3) { 
               llines(age., pred4$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., pred4$yhat.test.lower[1:K], col = 4, lwd = 2, lty = 3)
               llines(age., pred4$yhat.test.upper[1:K], col = 4, lwd = 2, lty = 3)
               ltext(3, 180, 'SRS: t,u', pos = 4)
           } else if(i == 4) { 
               llines(age., pred4$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., pred4$yhat.test.lower[K+1:K], col = 2, lwd = 2, lty = 3)
               llines(age., pred4$yhat.test.upper[K+1:K], col = 2, lwd = 2, lty = 3)
           } else if(i == 5) { 
               llines(age., pred0$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., pred0$yhat.test.lower[1:K], col = 4, lwd = 2, lty = 3)
               llines(age., pred0$yhat.test.upper[1:K], col = 4, lwd = 2, lty = 3)
               ltext(3, 180, 'MC: t,u', pos = 4)
           } else if(i == 6) { 
               llines(age., pred0$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., pred0$yhat.test.lower[K+1:K], col = 2, lwd = 2, lty = 3)
               llines(age., pred0$yhat.test.upper[K+1:K], col = 2, lwd = 2, lty = 3)
           } else if(i == 7) { 
               llines(age., pred2$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., pred2$yhat.test.lower[1:K], col = 4, lwd = 2, lty = 3)
               llines(age., pred2$yhat.test.upper[1:K], col = 4, lwd = 2, lty = 3)
               ltext(3, 180, 'SA: t,u,w(t,u)', pos = 4)
           } else if(i == 8) { 
               llines(age., pred2$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., pred2$yhat.test.lower[K+1:K], col = 2, lwd = 2, lty = 3)
               llines(age., pred2$yhat.test.upper[K+1:K], col = 2, lwd = 2, lty = 3)
           } else if(i == 9) { 
               llines(age., pred3$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., pred3$yhat.test.lower[1:K], col = 4, lwd = 2, lty = 3)
               llines(age., pred3$yhat.test.upper[1:K], col = 4, lwd = 2, lty = 3)
               ltext(3, 180, 'NN: t,u', pos = 4)
           } else if(i == 10) { 
               llines(age., pred3$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., pred3$yhat.test.lower[K+1:K], col = 2, lwd = 2, lty = 3)
               llines(age., pred3$yhat.test.upper[K+1:K], col = 2, lwd = 2, lty = 3)
           }
})
dev.off()
##dev.copy2pdf(file = 'FPD-final.pdf')

