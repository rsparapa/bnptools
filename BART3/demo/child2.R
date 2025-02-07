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

file. <- "child2NN5-pred.tu.RDS" 
if(!file.exists(file.)) {
    print(system.time(NN5.tu <- SHNN2(fit1, x.test.[ , 1:2], S = c('age', 'sex'))))
    saveRDS(NN5.tu, file.) 
} else {
    NN5.tu <- readRDS(file.) 
}
##     user   system  elapsed 
## 4862.598    2.152 4864.989 

NN5.tu. <- list()
NN5.tu.$yhat.test.mean <- apply(NN5.tu, 2, mean)
NN5.tu.$yhat.test.lower <- apply(NN5.tu, 2, quantile, probs = 0.025)
NN5.tu.$yhat.test.upper <- apply(NN5.tu, 2, quantile, probs = 0.975)

file. <- "child2NN5-pred.t.RDS" 
if(!file.exists(file.)) {
    print(system.time(NN5.t <- SHNN(fit1, age., S = 'age')))
    saveRDS(NN5.t, file.) 
} else {
    NN5.t <- readRDS(file.) 
}

NN5.t$yhat.test <- NN5.t$yhat.test+fit1$offset
NN5.t$yhat.test.mean <- NN5.t$yhat.test.mean+fit1$offset
NN5.t$yhat.test.lower <- NN5.t$yhat.test.lower+fit1$offset
NN5.t$yhat.test.upper <- NN5.t$yhat.test.upper+fit1$offset

file. <- "child2NN5-pred.w.RDS" 
if(!file.exists(file.)) {
    print(system.time(NN5.w <- SHNN(fit1, x.test.[ , 3], S = 'weight')))
    saveRDS(NN5.w, file.)
} else {
    NN5.w <- readRDS(file.) 
}

file. <- "child2NN5-pred.u.RDS" 
if(!file.exists(file.)) {
    print(system.time(NN5.u <- SHNN(fit1, 1:2, S = 'sex')))
    saveRDS(NN5.u, file.)
} else {
    NN5.u <- readRDS(file.) 
}

x.train <- bmxht[train, pick]
x.test  <- bmxht[!train, pick]

age. <- 2:17
K <- length(age.)
x.test. <- cbind(age., 0) 
x.test. <- rbind(x.test., x.test.)
x.test.[ , 2] <- rep(1:2, each = K)

(x.test. <- cbind(x.test.[ , 1:2], fit2$yhat.test.mean))


file. <- "child2SRS5-pred.tu.RDS" 
if(!file.exists(file.)) {
    print(system.time(SRS5.tu <- SHAP2(fit1, x.test.[ , 1:2], S = c('age', 'sex'))))
    saveRDS(SRS5.tu, file.) 
} else {
    SRS5.tu <- readRDS(file.) 
}
##     user   system  elapsed 
## 4862.598    2.152 4864.989 

SRS5.tu. <- list()
SRS5.tu.$yhat.test.mean <- apply(SRS5.tu, 2, mean)
SRS5.tu.$yhat.test.lower <- apply(SRS5.tu, 2, quantile, probs = 0.025)
SRS5.tu.$yhat.test.upper <- apply(SRS5.tu, 2, quantile, probs = 0.975)

file. <- "child2SRS5-pred.t.RDS" 
if(!file.exists(file.)) {
    print(system.time(SRS5.t <- SHAP(fit1, age., S = 'age')))
    saveRDS(SRS5.t, file.) 
} else {
    SRS5.t <- readRDS(file.) 
}

SRS5.t$yhat.test <- SRS5.t$yhat.test+fit1$offset
SRS5.t$yhat.test.mean <- SRS5.t$yhat.test.mean+fit1$offset
SRS5.t$yhat.test.lower <- SRS5.t$yhat.test.lower+fit1$offset
SRS5.t$yhat.test.upper <- SRS5.t$yhat.test.upper+fit1$offset

file. <- "child2SRS5-pred.w.RDS" 
if(!file.exists(file.)) {
    print(system.time(SRS5.w <- SHAP(fit1, x.test.[ , 3], S = 'weight')))
    saveRDS(SRS5.w, file.)
} else {
    SRS5.w <- readRDS(file.) 
}

file. <- "child2SRS5-pred.u.RDS" 
if(!file.exists(file.)) {
    print(system.time(SRS5.u <- SHAP(fit1, 1:2, S = 'sex')))
    saveRDS(SRS5.u, file.)
} else {
    SRS5.u <- readRDS(file.) 
}


df <- data.frame(sex = bmx$sex[train], age = bmx$age[train], height = bmx$height[train])
N <- nrow(df)
df <- rbind(df, df, df, df, df)
df$marginal <- rep(c('A', 'B', 'C', 'D', 'E'), each = N)
df$marginal <- factor(df$marginal)

pdf(file = 'SHAP-final5.pdf')
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
               llines(age., SRS5.u$yhat.test.mean[1]+SRS5.t$yhat.test.mean, col = 4, lwd = 2)
               llines(age., apply(SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test, 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test, 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
               ltext(1.5, 180, 'SRS t,u', pos = 4)
               ltext(15, 85, 'M', col = 4, pos = 4)
           } else if(i == 2) { 
               llines(age., SRS5.u$yhat.test.mean[2]+SRS5.t$yhat.test.mean, col = 2, lwd = 2)
               llines(age., apply(SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test, 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test, 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
               ltext(15, 85, 'F', col = 2, pos = 4)
           } else if(i == 3) { 
               ltext(1.5, 180, 'SRS t,u,t:u', pos = 4)
               llines(age., 2*SRS5.tu.$yhat.test.mean[1:K]+SRS5.u$yhat.test.mean[1]+SRS5.t$yhat.test.mean, col = 4, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , 1:K]+SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test, 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , 1:K]+SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test, 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 4) { 
               llines(age., 2*SRS5.tu.$yhat.test.mean[K+1:K]+SRS5.u$yhat.test.mean[2]+SRS5.t$yhat.test.mean, col = 2, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , K+1:K]+SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test, 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , K+1:K]+SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test, 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           } else if(i == 5) { 
               ltext(1.5, 180, 'NN t,u,t:u', pos = 4)
               llines(age., 2*NN5.tu.$yhat.test.mean[1:K]+NN5.u$yhat.test.mean[1]+NN5.t$yhat.test.mean, col = 4, lwd = 2)
               llines(age., apply(2*NN5.tu[ , 1:K]+NN5.u$yhat.test[ , 1]+NN5.t$yhat.test, 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*NN5.tu[ , 1:K]+NN5.u$yhat.test[ , 1]+NN5.t$yhat.test, 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 6) { 
               llines(age., 2*NN5.tu.$yhat.test.mean[K+1:K]+NN5.u$yhat.test.mean[2]+NN5.t$yhat.test.mean, col = 2, lwd = 2)
               llines(age., apply(2*NN5.tu[ , K+1:K]+NN5.u$yhat.test[ , 2]+NN5.t$yhat.test, 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*NN5.tu[ , K+1:K]+NN5.u$yhat.test[ , 2]+NN5.t$yhat.test, 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           } else if(i == 7) { 
               ltext(1.5, 180, 'SRS/SA t,u,t:u,w(t,u)', pos = 4)
               llines(age., 2*SRS5.tu.$yhat.test.mean[1:K]+SRS5.u$yhat.test.mean[1]+SRS5.t$yhat.test.mean+SRS5.w$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , 1:K]+SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test+SRS5.w$yhat.test[, 1:K], 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , 1:K]+SRS5.u$yhat.test[ , 1]+SRS5.t$yhat.test+SRS5.w$yhat.test[ , 1:K], 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 8) { 
               llines(age., 2*SRS5.tu.$yhat.test.mean[K+1:K]+SRS5.u$yhat.test.mean[2]+SRS5.t$yhat.test.mean+SRS5.w$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , K+1:K]+SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test+SRS5.w$yhat.test[, K+1:K], 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*SRS5.tu[ , K+1:K]+SRS5.u$yhat.test[ , 2]+SRS5.t$yhat.test+SRS5.w$yhat.test[ , K+1:K], 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           } else if(i == 9) { 
               ltext(1.5, 180, 'NN/SA t,u,t:u,w(t,u)', pos = 4)
               llines(age., 2*NN5.tu.$yhat.test.mean[1:K]+NN5.u$yhat.test.mean[1]+NN5.t$yhat.test.mean+NN5.w$yhat.test.mean[1:K], col = 4, lwd = 2)
               llines(age., apply(2*NN5.tu[ , 1:K]+NN5.u$yhat.test[ , 1]+NN5.t$yhat.test+NN5.w$yhat.test[, 1:K], 2, quantile, probs = 0.025), col = 4, lty = 3, lwd = 2)
               llines(age., apply(2*NN5.tu[ , 1:K]+NN5.u$yhat.test[ , 1]+NN5.t$yhat.test+NN5.w$yhat.test[ , 1:K], 2, quantile, probs = 0.975), col = 4, lty = 3, lwd = 2)
           } else if(i == 10) { 
               llines(age., 2*NN5.tu.$yhat.test.mean[K+1:K]+NN5.u$yhat.test.mean[2]+NN5.t$yhat.test.mean+NN5.w$yhat.test.mean[K+1:K], col = 2, lwd = 2)
               llines(age., apply(2*NN5.tu[ , K+1:K]+NN5.u$yhat.test[ , 2]+NN5.t$yhat.test+NN5.w$yhat.test[, K+1:K], 2, quantile, probs = 0.025), col = 2, lty = 3, lwd = 2)
               llines(age., apply(2*NN5.tu[ , K+1:K]+NN5.u$yhat.test[ , 2]+NN5.t$yhat.test+NN5.w$yhat.test[ , K+1:K], 2, quantile, probs = 0.975), col = 2, lty = 3, lwd = 2)
           }
})
dev.off()
     
