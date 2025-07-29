
library(BART3)
library(nftbart)

set.seed(34)
N <- 2500
L <- 130
Q <- N+L
P <- 20
x.train <- matrix(runif(P*N, -1, 1), nrow = N, ncol = P)
x.test  <- matrix(runif(P*Q, -1, 1), nrow = Q, ncol = P)
x <- seq(-1, 1.25, 0.25)
x.test[1:130, 1:5] <- 0
x.test[ 1:10, 1] <- x
x.test[11:20, 2] <- x
x.test[21:30, 3] <- x
h <- 30
for(i in 1:10)
    for(j in 1:10) {
        h <- h+1
        x.test[h, 4] <- x[i]
        x.test[h, 5] <- x[j]
}
print(h)

mu. <- function(Z) {
    Z <- rbind(Z)
    Z[ , 1]^3-Z[ , 2]^2+Z[ , 3]-Z[ , 4]*Z[ , 5]
}
sd. <- function(Z)  {
    Z <- rbind(Z)
    (Z[ , 6]^2)*(Z[ , 7]+1)
}
y.train <- exp(rnorm(N, mu.(x.train), sd.(x.train)))
c.train <- exp(rnorm(N, mu.(x.train), sd.(x.train)))
y.test  <- exp(rnorm(Q, mu.(x.test), sd.(x.test)))
summary(y.train)
summary(y.test[-(1:L)])
delta <- 1*(y.train<c.train)
times <- delta*y.train+(1-delta)*c.train
table(delta)/N

file. <- "ex7.rds"
K <- 25
XPtr <- TRUE
if(file.exists(file.)) {
    post <- readRDS(file.)
    XPtr <- FALSE
    if(mc.cores.openmp()>0) options(mc.cores = 8)
} else if(mc.cores.openmp()>0) {
    options(mc.cores = 8)
    set.seed(12)
    post <- nft(x.train, times = times, delta = delta,
                         x.test = x.test, K = K, nskip = 2000)
    saveRDS(post, file.)
    ## plot(post$sigma[ , 1], type = 'l', ylim = c(0, max(post$sigma)), 
    ##      ylab = expression(sigma))
    ## for(i in 2:8) lines(post$sigma[ , i], col = i)
    ## abline(v = 100, h = c(0, sd.), col = 8)
    ## abline(h = post$sigma.mean, lty = 2)
    ## check <- maxRhat(post$sigma., post$chains)
    ## acf(post$sigma., main = expression(sigma)) 
    ## points(check$rho)
    ## plot(post$accept[ , i], type = 'l', ylim = 0:1, ylab = 'MH')
    ## for(i in 2:8) lines(post$accept[ , i], col = i)
    ## abline(v = 100, h = 0:1, col = 8)
} else {
    set.seed(12)
    post <- nft(x.train, times = times, delta = delta, 
                x.test = x.test, K = K)
    saveRDS(post, file.)
    ## acf(post$sigma., main = expression(sigma)) 
}
  
B <- getOption('mc.cores', 1)
if(B>1) {
    s.max <- maxRhat(post$s.train.max, B)
    f.max <- maxRhat(post$f.train.max, B)
}

plot(post$s.train.max, type = 'l', 
     ylab = 'sd(x)', ylim = c(0, max(post$s.train.max)))
lines(post$s.train.min, col = 2)

acf(post$s.train.max)
if(B>1) points(s.max$rho)

plot(post$f.train.max, type = 'l', 
     ylab = 'mu(x)', 
     ylim = c(min(post$f.train.min), max(post$f.train.max)))
lines(post$f.train.min, col = 2)

acf(post$f.train.max)
if(B>1) points(f.max$rho)

print(c(test = Cindex(-post$f.test.mean[-(1:L)], y.test[-(1:L)]),
        train = Cindex(-post$f.train.mean, y.train),
        train. = Cindex(-post$f.train.mean, times, delta)))

## C <- apply(-post$f.train, 1, Cindex, y.train) ## very slow

plot(2:P, cumsum(post$f.varprob)[-1], 
     xlab = 'x', ylab = 's for mu(x)', type = 'p', 
     xlim = c(1, P), ylim = 0:1, col = 4) 
points(1:P, post$f.varprob, col = 1+(post$f.varprob>(1/P)))
abline(h = c(0, 1/P, 1), col = c(1, 8, 1), lty = c(1, 2, 1))
text(P, 1/P, '1/P', pos = 3)
abline(v = 5, col = 8)
legend('right', legend = c('Active', 'Inactive', 'Cumulative'),
       col = c(2, 1, 4), pch = 1)

plot(2:P, cumsum(post$s.varprob)[-1], 
     xlab = 'x', ylab = 's for sigma(x)', type = 'p', 
     xlim = c(1, P), ylim = 0:1, col = 4) 
points(1:P, post$s.varprob, col = 1+(post$s.varprob>(1/P)))
abline(h = c(0, 1/P, 1), col = c(1, 8, 1), lty = c(1, 2, 1))
text(P, 1/P, '1/P', pos = 3)
abline(v = c(5, 7), col = 8)
legend('right', legend = c('Active', 'Inactive', 'Cumulative'),
       col = c(2, 1, 4), pch = 1)

times. <- exp(seq(quantile(log(y.train), 0.025),
                  quantile(log(y.train), 0.975), length.out = K))
summary(times.)

pred <- NULL
## this takes a while
## pred <- predict(post, x.test, XPtr = XPtr, events = times., K = K, na.rm = TRUE)

j <- 1
## j <- 2
## j <- 3
par(mfrow = c(3, 3))
for(i in 1:9) {
    h <- (i-1)*K+1:K
    k <- (j-1)*10+i
    plot(c(0, times.), 
         c(1, pnorm(log(times.), mu.(x.test[k, ]), sd.(x.test[k, ]), FALSE)), 
         type = 'l', ylim = 0:1, xlab = 't', ylab = 'S(t|x)')
    if(length(pred) == 0) {
    lines(c(0, times.), c(1, pnorm(log(times.), 
          post$f.test.mean[k], post$s.test.mean[k], FALSE)), col = 2)
    } else {
        lines(c(0, times.), c(1, pred$surv.test.mean[h]), col = 2)
        lines(c(0, times.), c(1, pred$surv.test.lower[h]), col = 4)
        lines(c(0, times.), c(1, pred$surv.test.upper[h]), col = 4)
    }
    abline(v = 0, h = 0:1, col = 8)
    text(3.5, 0.75, c(expression(x[1]), expression(x[2]), 
                    expression(x[3]))[j], pos = 2)
    text(3.5, 0.75, '=')
    text(3.5, 0.75, x.test[k, j], pos = 4)
}
par(mfrow = c(1, 1))

