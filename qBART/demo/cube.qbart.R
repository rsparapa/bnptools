
library(qBART)

set.seed(33120)
#simulate data
p <- 0.7  #not cured
n <- 2000  #total subj
status <- sample(0:1, n, replace = TRUE, prob = c((1-p),p))
x1 <- rnorm(n, mean = 45, sd = 10)
x2 <- sample(0:1, n, replace = TRUE)
b0 <- 1; b1 <- 0.02; b2 <- 0.13
xb <- b0 + b1*x1 + b2*x2
ltime <- rnorm(n, mean = xb)  #log(real_time)~Normal
time <- exp(ltime)
ctime <- rexp(n, rate = 0.01)  #censoring_time~Exp(0.01)
simcure <- data.frame(age = x1, gender = x2, c_time = ctime, r_time = time)
simcure <- simcure[order(simcure$r_time, decreasing = TRUE),]  #ordered real time decreasing
simcure$Ncure <- sort(status)  #larger real_time means cured
simcure[simcure$c_time > 30,]$c_time <- 30  #admin censoring
simcure$obstime <- ifelse(simcure$Ncure == 1, ifelse(simcure$c_time < simcure$r_time, simcure$c_time, simcure$r_time), simcure$c_time)
simcure$event <- ifelse(simcure$obstime == simcure$c_time, 0, 1)

x.train <- simcure[,c("age", "gender")]
times <- simcure$obstime
delta <- simcure$event

post1 <- qbart(x.train1=x.train, x.train2=x.train, times=times, delta=delta, x.test1=x.train)

library(flexsurvcure)

bcx <- bc$group
bct <- bc$rectime
bcd <- bc$censrec
postbc <- qbart(x.train1=bcx, x.train2=bcx, times=bct, delta=bcd)
str(postbc)
mean(postbc$prob.train)

str(t(postbc$x.train1))
temp <- cbind(bc,t(postbc$x.train1))
mix1 <- flexsurvcure(Surv(rectime,censrec)~ X2 + X3 + meanlog(X2) + meanlog(X3), data=temp, dist="lnorm", link="logistic", mixture = T)

mixture = flexsurvcure(Surv(rectime,censrec)~ group + meanlog(group), data=bc, dist="lnorm", link="logistic", mixture = T)
print(mixture)

str(lung)
lx <- lung[,-(1:3)]
lt <- lung$time
ld <- lung$status-1
postl <- qbart(x.train1=lx, x.train2=lx, times=lt, delta=ld)

plot(postl$times, postl$surv.train.mean[1:K], type='s')

par(mfrow=c(1,2))
kmfit <- survfit(Surv(time, status)~1, data=lung)
plot(kmfit)

K <- postbc$K
N <- nrow(lung)
for(j in 1:K) {
    postl$surv.test.025 <- matrix(nrow=ndpost, ncol=K)
    postl$surv.test.025[ , j] <-
                apply(postl$surv.train[ , seq(j, K*N, K)],
                      1, mean)
    postl$surv.test.025 <- apply(postl$surv.test.025,2,mean)
    plot(c(0, postl$times),
              c(1, postl$surv.test.025),
              col=2, type='s', ylim=0:1, xlab='t', ylab='S(t, x)')
}
