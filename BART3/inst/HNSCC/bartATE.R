
data(hnscc, package = 'BART3')

pred <- readRDS("bart.rds")

K <- pred$K
N <- length(hnscc$DFS)

## the 24 month time point is exceeded in interval 66 (24.04731)
pred$times[65]
h <- seq(65, N*K, K)
surv.est <- pred$surv.test[ , c(h, N*K+h)]

## start with survival
ITE1 <- surv.est[, 1:N]
ITE0 <- surv.est[, -(1:N)]

observed <- cbind(hnscc$Study == 1, hnscc$Study == 0)

## account for time=24 of observed potential outcomes
failure24 <- (hnscc$DFS_Time <= 24 & hnscc$DFS == 1)
table(failure24, useNA = 'ifany')

success24 <- (hnscc$DFS_Time >= 24)
table(success24, useNA = 'ifany')

table(failure24 & success24, useNA = 'ifany') ## you cannot be both

ITE1[observed[ , 1] & failure24, 1] <- 0
ITE1[observed[ , 1] & success24, 1] <- 1

ITE0[observed[ , 2] & failure24, 2] <- 0
ITE0[observed[ , 2] & success24, 2] <- 1

ATE <- apply(ITE1-ITE0, 1, mean)

quantile(ATE, probs = 0.025)
mean(ATE)
quantile(ATE, probs = 0.975)

