
library(BART3)
library(hbart)

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
sd. <- function(Z) (Z[ , 6]^2)*(Z[ , 7]+1)
y.train <- rnorm(N, mu.(x.train), sd.(x.train))
y.test  <- rnorm(Q, mu.(x.test), sd.(x.test))
summary(y.train)
summary(y.test[-(1:L)])

XPtr <- TRUE
file. <- "ex5.rds"
if(file.exists(file.)) {
    post <- readRDS(file.)
    XPtr <- FALSE
    if(mc.cores.openmp()>0) options(mc.cores = 8)
} else {
    if(mc.cores.openmp()>0) {
        options(mc.cores = 8)
    } else {
        options(mc.cores = 1)
    }
    set.seed(12)
    post <- hbart(x.train, y.train, x.test)
    saveRDS(post, file.)
}
   
plot(post$f.train.mean, y.train, asp = TRUE,
     xlab = 'f(x)', ylab = 'y')
text(-3, 2, expression(R^2), pos = 2)
text(-3, 2, '=')
text(-3, 2, round(
              cor(y.train, post$f.train.mean)^2,
              digits = 3), pos = 4)
abline(a = 0, b = 1, col = 8)

plot(post$f.test.mean[-(1:L)], y.test[-(1:L)], asp = TRUE, 
     xlab = 'f(x)', ylab = 'y')
text(-3, 2, expression(R^2), pos = 2)
text(-3, 2, '=')
text(-3, 2, round(
              cor(y.test[-(1:L)], post$f.test.mean[-(1:L)])^2,
              digits = 3), pos = 4)
abline(a = 0, b = 1, col = 8)

plot(post$s.train.mean, sd.(x.train), asp = TRUE,
     xlab = 'HBART', ylab = 's(x)')
text(0, 1.5, expression(R^2), pos = 2)
text(0, 1.5, '=')
text(0, 1.5, round(
              cor(sd.(x.train), post$s.train.mean)^2,
              digits = 3), pos = 4)
abline(a = 0, b = 1, col = 8)

plot(post$f.test.mean[-(1:L)], y.test[-(1:L)], asp = TRUE, 
     xlab = 'f(x)', ylab = 'y')
text(-3, 2, expression(R^2), pos = 2)
text(-3, 2, '=')
text(-3, 2, round(
              cor(y.test[-(1:L)], post$f.test.mean[-(1:L)])^2,
              digits = 3), pos = 4)
abline(a = 0, b = 1, col = 8)

plot(post$s.test.mean[-(1:L)], sd.(x.test)[-(1:L)], asp = TRUE,
     xlab = 'HBART', ylab = 's(x)')
text(0, 1.5, expression(R^2), pos = 2)
text(0, 1.5, '=')
text(0, 1.5, round(
              cor(sd.(x.test)[-(1:L)], post$s.test.mean[-(1:L)])^2,
              digits = 3), pos = 4)
abline(a = 0, b = 1, col = 8)

## plot(2:P, cumsum(post$varprob.mean)[-1], 
##      xlab = 'x', ylab = 's', type = 'p', 
##      xlim = c(1, P), ylim = 0:1, col = 4) 
## points(1:P, post$varprob.mean, col = 1+(post$varprob.mean>(1/P)))
## abline(h = c(0, 1/P, 1), col = c(1, 8, 1), lty = c(1, 2, 1))
## text(P, 1/P, '1/P', pos = 3)
## abline(v = 5, col = 8)
## legend('right', legend = c('Active', 'Inactive', 'Cumulative'),
##        col = c(2, 1, 4), pch = 1)

plot(x.test[1:10, 1], mu.(x.test[1:10, ]), type = 'n', 
     ylim = c(min(post$yhat.test.lower[1:L]),
              max(post$yhat.test.upper[1:L])),
     xlab = 'x', ylab = 'f(x)')
abline(v = 1, col = 8)
legend('topleft', legend = c('+cubic', '-quadratic', '+linear'),
       lwd = 2, col = c(1, 2, 4))
for(i in 1:3) {
    h <- 2^(i-1)
    j <- (1+(i-1)*10):(i*10)
    X <- matrix(0, nrow = L, ncol = P)
    X[ , i] <- seq(-1, 1.25, length.out = L)
    lines(X[ , i], mu.(X), lwd=2, col=h)
    lines(x.test[j, i], post$yhat.test.upper[j], lty=2, col=h)
    lines(x.test[j, i], post$yhat.test.lower[j], lty=2, col=h)
    points(x.test[j[9:10], i], post$yhat.test.mean[j[9:10]], col=h)
}
 
print(cor(y.test[1:L], post$yhat.test.mean[1:L])^2)

levels. <- quantile(-outer(x[-10], x[-10]), (1:4)/5)
contour(x, x, -outer(x, x), levels = levels.)
abline(v = 1, h = 1, col = 8)
z <- matrix(nrow = 10, ncol = 10)
h <- 30
for(i in 1:10)
    for(j in 1:10) {
        h <- h+1
        z[i, j] <- post$yhat.test.mean[h] ## x4:i, x5:j
}
contour(x, x, z, add = TRUE, col = 2, levels = levels.)
