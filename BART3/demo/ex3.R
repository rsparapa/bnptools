
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

mu. <- function(Z) Z[ , 1]^3-Z[ , 2]^2+Z[ , 3]-Z[ , 4]*Z[ , 5]
sd. <- 0.25
y.train <- rnorm(N, mu.(x.train), sd.)
y.test  <- rnorm(Q, mu.(x.test), sd.)
summary(y.train)
summary(y.test[-(1:L)])

Y.train <- rep(0, N)
Y.train[y.train > 0] <- 3
Y.train[y.train <= 0 & y.train> -1] <- 2
Y.train[y.train <= -1] <- 1
table(Y.train)/N

Y.test <- rep(0, Q)
Y.test[y.test > 0] <- 3
Y.test[y.test <= 0 & y.test> -1] <- 2
Y.test[y.test <= -1] <- 1
table(Y.test)/Q

prob.train <- matrix(nrow = N, ncol = 3)
prob.train[ , 1] <- pnorm(-1, mu.(x.train), sd.)
prob.train[ , 2] <- pnorm( 0, mu.(x.train), sd.)-prob.train[ , 1]
prob.train[ , 3] <- pnorm( 0, mu.(x.train), sd., FALSE)
summary(prob.train[ , 1]+prob.train[ , 2]+prob.train[ , 3])

prob.test <- matrix(nrow = Q, ncol = 3)
prob.test[ , 1] <- pnorm(-1, mu.(x.test), sd.)
prob.test[ , 2] <- pnorm( 0, mu.(x.test), sd.)-prob.test[ , 1]
prob.test[ , 3] <- pnorm( 0, mu.(x.test), sd., FALSE)
summary(prob.test[ , 1]+prob.test[ , 2]+prob.test[ , 3])

post <- list()
file. <- "ex3.rds"
if(file.exists(file.)) {
    post <- readRDS(file.)
} else if(.Platform$OS.type == 'unix') {
    options(mc.cores = 8)
    for(i in 1:3)
        post[[i]] <- mc.gbart(x.train, 1*(Y.train == i), x.test, 
                         type = 'pbart', sparse = TRUE, seed = 21)
    saveRDS(post, file.)
} else {
    set.seed(12)
    for(i in 1:3)
        post[[i]] <- gbart(x.train, 1*(Y.train == i), x.test, 
                         type = 'pbart', sparse = TRUE)
    saveRDS(post, file.)
}
   
prob.test.mean <- cbind(post[[1]]$prob.test.mean,
                        post[[2]]$prob.test.mean,
                        post[[3]]$prob.test.mean)

Yhat.test <- t(apply(prob.test.mean, 1, order))[ , 3]
table(Y.test[-(1:L)])
table(Y.test[-(1:L)], Yhat.test[-(1:L)])/(Q-L)

for(i in 1:3) {
    print(c(auROC = 
          mean(outer(post[[i]]$prob.test.mean[-(1:L)][Y.test[-(1:L)] != i], 
                     post[[i]]$prob.test.mean[-(1:L)][Y.test[-(1:L)] == i], 
                     '<'))))
    post[[i]]$yhat.train.max <- apply(post[[i]]$yhat.train, 1, max)
    post[[i]]$yhat.train.min <- apply(post[[i]]$yhat.train, 1, min)
    if(i == 1) {
        plot(post[[i]]$yhat.train.max, type = 'l', ylab = 'max',
             ylim = range(c(post[[1]]$yhat.train, post[[2]]$yhat.train, 
                            post[[3]]$yhat.train)))
    } else {
        lines(post[[i]]$yhat.train.max, col = 2^(i-1))
    }
    lines(post[[i]]$yhat.train.min, col = 2^(i-1), lty = 2)
}

par(mfrow = c(3, 2))
for(i in 1:3) {
plot(post[[i]]$prob.test.mean[-(1:L)], prob.test[-(1:L), i], 
     asp = TRUE, ylab = paste0('P(Y=', i, ')'), cex = 0.5,
     xlab = 'BART', xlim = 0:1, ylim = 0:1)
abline(a = 0, b = 1, col = 8, v = 0:1, h = 0:1)
text(-0.5, 0.9, expression(R^2), pos = 2)
text(-0.5, 0.9, '=')
text(-0.5, 0.9, round(
              cor(prob.test[-(1:L), i], 
                  post[[i]]$prob.test.mean[-(1:L)])^2,
              digits = 3), pos = 4)
plot(1:P, post[[i]]$varprob.mean, xlab = 'x', ylab = 's', type = 'p', 
ylim = 0:1, col = 1+(post[[i]]$varprob.mean>(1/P)))
points(2:P, cumsum(post[[i]]$varprob.mean)[-1], col = 4)
abline(h = c(0, 1/P, 1), col = c(1, 8, 1), lty = c(1, 2, 1))
abline(v = 5, col = 8)
text(P, 1/P, '1/P', pos = 3)
}
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
for(i in 1:3) {
    print(cor(prob.test[1:L, i], post[[i]]$prob.test.mean[1:L])^2)
    plot(x.test[1:10, 1], prob.test[1:10, i], type = 'n', 
         ylim = 0:1, xlab = 'x', ylab = paste0('P(Y=', i, ')'))
    abline(v = 1, col = 8)
    if(i == 1)
        legend('top', legend = c('+cubic', '-quadratic', '+linear'),
           lwd = 2, col = c(1, 2, 4))
    for(k in 1:3) {
        h <- 2^(k-1)
        j <- (1+(k-1)*10):(k*10)
        X <- matrix(0, nrow = L, ncol = P)
        X[ , k] <- seq(-1, 1.25, length.out = L)
        prob.test. <- double(L)
        prob.test. <- pnorm(-1, mu.(X), sd.)
        if(i == 2) prob.test. <- pnorm( 0, mu.(X), sd.)-prob.test.
        if(i == 3) prob.test. <- pnorm( 0, mu.(X), sd., FALSE)
        lines(X[ , k], prob.test., lwd=2, col=h)
        lines(x.test[j, k], post[[i]]$prob.test.upper[j], lty=2, col=h)
        lines(x.test[j, k], post[[i]]$prob.test.lower[j], lty=2, col=h)
        points(x.test[j[9:10], k], post[[i]]$prob.test.mean[j[9:10]], col=h)
    }
}
par(mfrow = c(1, 1))
 

