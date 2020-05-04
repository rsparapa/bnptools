library(qBART)

set.seed(33120)
#simulate data
p <- 0.7  #not cured
n <- 2000  #total subj
status <- sample(0:1, n, replace = TRUE, prob = c((1-p),p))
x1 <- rnorm(n, mean = 0, sd = 4)
#x2 <- sample(0:1, n, replace = TRUE) gender = x2,
b0 <- 1; b1 <- 0.5; #b2 <- 0.13
xb <- b0 + b1*x1 #+ b2*x2
ltime <- rnorm(n, mean = xb)  #log(real_time)~Normal
time <- exp(ltime)
ctime <- rexp(n, rate = 0.01)  #censoring_time~Exp(0.01)
simcure <- data.frame(age = x1,  c_time = ctime, r_time = time)
#simcure <- simcure[order(simcure$r_time, decreasing = TRUE),]  #ordered real time decreasing
simcure$Ncure <- status  #larger real_time means cured
cuttime <- quantile(time, probs=0.8) 
simcure[simcure$c_time > cuttime,]$c_time <- cuttime  #admin censoring
simcure$r_time <- ifelse(simcure$Ncure == 1, simcure$r_time, Inf)
simcure$obstime <- ifelse(simcure$c_time < simcure$r_time, simcure$c_time, simcure$r_time)
simcure$event <- ifelse(simcure$obstime == simcure$c_time, 0, 1)

x.train <- simcure$age
times <- simcure$obstime
delta <- simcure$event

notcure <- simcure$Ncure==1
tdata <- simcure[notcure,]
library(BART3)
post <- abart(x.train=tdata$age, times=tdata$obstime, delta=tdata$event)
#plot(xb[notcure], post$yhat.train.mean)
plot(xb[notcure], ltime[notcure])
par(mfrow=c(1,2))
plot(ltime[notcure],post$yhat.train.mean,pch=20,main="abart")
abline(a=0,b=1,col=2)

post1 <- qbart(x.train1=x.train, x.train2=x.train, times=times, delta=delta)
#str(post1)
plot(ltime[notcure],post1$y2hat.train.mean[notcure], pch=20,main="qbart")
abline(a=0,b=1,col=2)
mean(post1$prob.train)

K <- post1$K
total.surv.train <- matrix(0,nrow=n,ncol=K)
for (i in 1:n){
    total.surv.train[i,] <- post1$surv.train.mean[(i-1)*K+(1:K)]
}
total.surv.mean <- apply(total.surv.train, 2, mean)
## par(mfrow=c(1,2))
plot(c(0,post1$times),c(1,total.surv.mean),ylim=c(0,1),type='s')
kmfit1 <- survfit(Surv(times, delta)~1)
plot(kmfit1)

library(flexsurvcure)

bcx <- bc$group
bct <- bc$rectime
bcd <- bc$censrec
postbc <- qbart(x.train1=bcx, x.train2=bcx, times=bct, delta=bcd)
## str(postbc)
mean(postbc$prob.train)
n <- nrow(bc)
K <- postbc$K
total.surv.train <- matrix(0,nrow=n,ncol=K)
for (i in 1:n){
    total.surv.train[i,] <- postbc$surv.train.mean[(i-1)*K+(1:K)]
}
total.surv.mean <- apply(total.surv.train, 2, mean)
par(mfrow=c(1,2))
plot(c(0,postbc$times),c(1,total.surv.mean),ylim=c(0,1), type='s')
sfit <- survfit(Surv(rectime, censrec)~1, data=bc)
plot(sfit)

mixture = flexsurvcure(Surv(rectime,censrec)~ group + meanlog(group), data=bc, dist="lnorm", link="logistic", mixture = T)
print(mixture)

str(lung)
lx <- lung[,-(1:3)]
lt <- lung$time
ld <- lung$status-1
postl <- qbart(x.train1=lx, x.train2=lx, times=lt, delta=ld)
str(postl)
mean(postl$prob.train)
mix <- flexsurvcure(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + meal.cal + wt.loss, data=lung, dist="lnorm")

n <- nrow(lung)
K <- postl$K
total.surv.train <- matrix(0,nrow=n,ncol=K)
for (i in 1:n){
    total.surv.train[i,] <- postl$surv.train.mean[(i-1)*K+(1:K)]
}
total.surv.mean <- apply(total.surv.train, 2, mean)
par(mfrow=c(1,2))
plot(c(0,postl$times),c(1,total.surv.mean),ylim=c(0,1),type='s')
kmfit <- survfit(Surv(time, status)~1, data=lung)
plot(kmfit)

