
library(BART3)

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

sd. <- 0.5
y.train <- exp(rnorm(N, mu.(x.train), sd.))
c.train <- exp(rnorm(N, mu.(x.train), sd.))
y.test  <- exp(rnorm(Q, mu.(x.test), sd.))
summary(y.train)
summary(y.test[-(1:L)])
delta <- 1*(y.train<c.train)
times <- delta*y.train+(1-delta)*c.train
table(delta)/N

file. <- "ex6.rds"
K <- 50
if(file.exists(file.)) {
    post <- readRDS(file.)
} else if(.Platform$OS.type == 'unix') {
    options(mc.cores = 8)
    post <- mc.surv.bart(x.train, times = times, delta = delta,
                         x.test = x.test, K = K, sparse = TRUE, seed = 21)
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
    post <- surv.bart(x.train, times = times, delta = delta, x.test = x.test, 
                      K = K, sparse = TRUE)
    saveRDS(post, file.)
    ## acf(post$sigma., main = expression(sigma)) 
}
  
C <- 0
for(i in 1:K) 
    C[i] <- Cindex(1-post$surv.test.mean[L*K+seq(i, N*K, K)], y.test[-(1:L)])
summary(C)

plot(1:P, cumsum(post$varprob.mean)[-1], 
     xlab = 'x', ylab = 's', type = 'p', 
     xlim = c(0, P), ylim = 0:1, col = 4) 
points(0:P, post$varprob.mean, col = 1+(post$varprob.mean>(1/(P+1))))
abline(h = c(0, 1/(P+1), 1), col = c(1, 8, 1), lty = c(1, 2, 1))
text(0, post$varprob.mean[1], 'Time', pos = 4)
text(P, 1/(P+1), '1/P', pos = 3)
abline(v = 5, col = 8)
legend('right', legend = c('Active', 'Inactive', 'Cumulative'),
       col = c(2, 1, 4), pch = 1)

times. <- seq(min(post$times), max(post$times), length.out = L)
j <- 1
## j <- 2
## j <- 3
par(mfrow = c(3, 3))
for(i in 1:9) {
    h <- (i-1)*K+1:K
    k <- (j-1)*10+i
    plot(c(0, times.), 
         c(1, pnorm(log(times.), mu.(x.test[k, ]), sd., FALSE)), 
         type = 'l', ylim = 0:1, xlab = 't', ylab = 'S(t|x)')
    ##lines(c(0, post$times), c(1, post$surv.test.mean[h]), col = 2)
    lines(c(0, post$times), c(1, post$surv.test.lower[h]), col = 2)
    lines(c(0, post$times), c(1, post$surv.test.upper[h]), col = 2)
    abline(v = 0, h = 0:1, col = 8)
    text(5, 0.75, c(expression(x[1]), expression(x[2]), 
                    expression(x[3]))[j], pos = 2)
    text(5, 0.75, '=')
    text(5, 0.75, x.test[k, j], pos = 4)
}
par(mfrow = c(1, 1))

