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

acf(fit1$sigma., ci.type = 'ma')

C <- ncol(fit1$sigma)
if(C>1) { ## requires parallel chains typically multi-threaded

    check <- maxRhat(fit1$sigma., 8)
    print(check$maxRhat)

    plot(check$splitrho[1:30], ylim = c(-0.1, 1), type = 'h',
         ylab = 'Auto-correlation', sub = 'Split Rhat rho')
    abline(h = 0:1, col = 8)

    acf(fit1$sigma., ci.type = 'ma', col = 2, lwd = 2, main = '')
    points(check$splitrho[1:30])
    ##dev.copy2pdf(file='growth0-acf.pdf')
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

}

pred0 <- predict(fit1, x.test)
pred0.mean <- apply(pred0, 2, mean)
print(c(test.Rsqr = cor(pred0.mean, fit1$yhat.test.mean)^2))
plot(pred0.mean, fit1$yhat.test.mean)

library(lattice)
K <- 30
(ndpost <- fit1$ndpost)
df <- data.frame(pred = c(pred0[ , 1:K]), 
                 SEQN = rep(bmx$SEQN[!train][1:K], each = ndpost))
densityplot(~pred|SEQN, data = df, layout = c(6, 5), plot.points = FALSE,
            ##strip = strip.custom(strip.names = FALSE), as.table = TRUE)
            strip = FALSE, as.table = TRUE,
            panel=function(...) {
                i <- panel.number()
                panel.abline(v = bmx$BMXHT[!train][i], col = 2)
                panel.abline(h = 0, col = 8)
                panel.densityplot(...)
            })
