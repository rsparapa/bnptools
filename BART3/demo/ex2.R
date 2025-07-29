
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

Y.train <- 1*(y.train > 0)
table(Y.train)/N

Y.test <- 1*(y.test > 0)
table(Y.test)/Q

prob.train <- pnorm(0, mu.(x.train), sd., FALSE)
prob.test  <- pnorm(0, mu.(x.test), sd., FALSE)

file. <- "ex2.rds"
if(file.exists(file.)) {
    post <- readRDS(file.)
} else if(.Platform$OS.type == 'unix') {
    options(mc.cores = 8)
    post <- mc.gbart(x.train, Y.train, x.test, 
                     type = 'pbart', sparse = TRUE, seed = 21)
    saveRDS(post, file.)
} else {
    set.seed(12)
    post <- gbart(x.train, Y.train, x.test, 
                  type = 'pbart', sparse = TRUE)
    saveRDS(post, file.)
}
 
yhat.train.max <- apply(post$yhat.train, 1, max)
yhat.train.min <- apply(post$yhat.train, 1, min)
plot(yhat.train.max, type = 'l', ylab = 'Max/Min',
     ylim = range(post$yhat.train))
lines(yhat.train.min, col = 2)

if(.Platform$OS.type == 'unix') {
    check.max <- maxRhat(yhat.train.max, C = 8)
    acf(yhat.train.max, main = 'Max')
    points(check.max$rho)
    readline()
    check.min <- maxRhat(yhat.train.min, C = 8)
    acf(yhat.train.min, main = 'Min')
    points(check.min$rho)
} else {
    acf(yhat.train.max, main = 'Max')
    readline()
    acf(yhat.train.min, main = 'Min')
}

prob.test.mean <- post$prob.test.mean[-(1:L)]
p <- seq(0, 1, length.out = L)
Y.test. <- Y.test[-(1:L)]
sens <- 0
spec <- 0
for(i in 1:L) {
    a <- sum((Y.test. == 0) & (prob.test.mean<p[i]))
    b <- sum((Y.test. == 1) & (prob.test.mean<p[i]))
    c <- sum((Y.test. == 0) & (prob.test.mean >= p[i]))
    d <- sum((Y.test. == 1) & (prob.test.mean >= p[i]))
    sens[i] <- d/(d+b)
    spec[i] <- a/(a+c)
}
plot(1-spec, sens, type = 'l',
     xlab = '1-Specificity', ylab = 'Sensitivity')
abline(a = 0, b = 1)
text(0.25, 0.75, paste0('Area under ROC = ', round(
          mean(outer(prob.test.mean[Y.test. == 0], 
                     prob.test.mean[Y.test. == 1], 
                     '<')), 3))) 

plot(prob.test.mean, prob.test[-(1:L)], 
     asp = TRUE, ylab = 'P(Y=1)', cex = 0.5,
     xlab = 'BART', xlim = 0:1, ylim = 0:1)
abline(a = 0, b = 1, col = 8, v = 0:1, h = 0:1)
text(0, 0.9, expression(R^2), pos = 2)
text(0, 0.9, '=')
text(0, 0.9, round(cor(prob.test[-(1:L)], prob.test.mean)^2,
              digits = 3), pos = 4)

plot(1:P, post$varprob.mean, xlab = 'x', ylab = 's', type = 'p', 
     ylim = 0:1, col = 1+(post$varprob.mean>(1/P)))
points(2:P, cumsum(post$varprob.mean)[-1], col = 4)
abline(h = c(0, 1/P, 1), col = c(1, 8, 1), lty = c(1, 2, 1))
abline(v = 5, col = 8)
text(P, 1/P, '1/P', pos = 3)

plot(x.test[1:10, 1], prob.test[1:10], type = 'n', 
     ylim = 0:1, xlab = 'x', ylab = 'P(Y=1)')
abline(v = 1, col = 8)
legend('topleft', legend = c('+cubic', '-quadratic', '+linear'),
       lwd = 2, col = c(1, 2, 4))
for(k in 1:3) {
    h <- 2^(k-1)
    j <- (1+(k-1)*10):(k*10)
    X <- matrix(0, nrow = L, ncol = P)
    X[ , k] <- seq(-1, 1.25, length.out = L)
    lines(X[ , k], pnorm(0, mu.(X), sd., FALSE), lwd=2, col=h)
    lines(x.test[j, k], post$prob.test.upper[j], lty=2, col=h)
    lines(x.test[j, k], post$prob.test.lower[j], lty=2, col=h)
    points(x.test[j[9:10], k], post$prob.test.mean[j[9:10]], col=h)
}
 
levels. <- quantile(outer(x[-10], x[-10], 
           function(a, b) pnorm(0, -a*b, sd., FALSE)), (1:4)/5)
contour(x, x, outer(x, x, 
           function(a, b) pnorm(0, -a*b, sd., FALSE)), levels = levels.)
abline(v = 1, h = 1, col = 8)
z <- matrix(nrow = 10, ncol = 10)
h <- 30
for(i in 1:10)
    for(j in 1:10) {
        h <- h+1
        z[i, j] <- post$prob.test.mean[h] ## x4:i, x5:j
}
contour(x, x, z, add = TRUE, col = 2, levels = levels.)
