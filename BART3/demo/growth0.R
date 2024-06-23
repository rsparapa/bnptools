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

pred0 <- predict(fit1, x.test)
pred0.mean <- apply(pred0, 2, mean)
print(c(test.Rsqr = cor(pred0.mean, fit1$yhat.test.mean)^2))
pred0. <- density(pred0[ , 1])
plot(pred0., main='First child: SEQN=1', xlab='Height (cm)')
abline(h = 0, v = quantile(pred0[ , 1], probs = c(0.025, 0.975)), col = 8)

C <- ncol(fit1$sigma)
if(C>1) { ## requires parallel chains typically multi-thread
    check <- maxRhat(fit1$sigma., 8)
    print(check$maxRhat)
    sigma.mean <- mean(fit1$sigma.) 
    sigma.max <- max(c(fit1$sigma))
    sigma.y <- mean(c(sigma.max, sigma.mean))
    plot(fit1$sigma[ , 1], type = 'l', 
         ylab = expression(sigma))
    for(i in 2:C) lines(fit1$sigma[ , i], col = i)
    abline(v = c(0, 100), h = sigma.mean, col = 8)
    text(50, sigma.y, "Burn-in discarded")
    text(160, sigma.y, "Converged posterior samples kept")
    ##dev.copy2pdf(file='growth0-sigma.pdf')
    plot(check$splitrho, ylim = c(-1, 1), type = 'h',
         ylab = 'Auto-correlation', sub = 'Split Rhat rho')
    abline(h = (-1):1, col = 8)
    ##dev.copy2pdf(file='growth0-sigma.pdf')
}

