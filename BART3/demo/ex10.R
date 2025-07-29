
library(BART3)
library(LearnBayes)

set.seed(34)
N <- 2500
L <- 130
Q <- N+L
P <- 20
sd. <- 0.25
rho. <- 0.7
Sigma <- matrix(sd.^2, nrow = P, ncol = P)-diag(sd.^2, nrow = P, ncol = P)
Sigma <- Sigma*rho.+diag(sd.^2, nrow = P, ncol = P)
X.train <- rmnorm(N, rep(0, P), Sigma)
X.train[X.train< -1] <- -1
X.train[X.train>1] <- 1
NA.train <- matrix(0, nrow = N, ncol = P)
for(i in 1:5) {
    NA.train[ , 6-i] <- rbinom(N, 1, i/10) 
    print(table(NA.train[ , 6-i])/N)
}
x.train <- X.train
x.train[which(NA.train == 1)] <- NA

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
sd. <- 0.5
y.train <- rnorm(N, mu.(X.train), sd.)
y.test  <- rnorm(Q, mu.(x.test), sd.)
summary(y.train)
summary(y.test[-(1:L)])

saveRDS(x.train, "ex10-x.train.rds")
saveRDS(y.train, "ex10-y.train.rds")

## file. <- "ex10.rds"
## if(file.exists(file.)) {
##     post <- readRDS(file.)
## } else {
##     set.seed(12)
##     post <- BartMI(paste0(getwd(), '/'), xx = x.train, yy = y.train, 
##                    datatype = rep(0, P), type = 1)
## }
   
