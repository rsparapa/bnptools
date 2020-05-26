library(qBART)

set.seed(33120)
#simulate data; p independent
p <- 0.5  #not cured
n <- 2000  #total subj
status <- sample(0:1, n, replace = TRUE, prob = c((1-p),p))
x1 <- rnorm(n, mean = 0, sd = 4)
#x2 <- sample(0:1, n, replace = TRUE)
b0 <- 1; b1 <- 0.5; #b2 <- 0.13
xb <- b0 + b1*x1 #+ b2*x2
ltime <- rnorm(n, mean = xb)  #log(real_time)~Normal
time <- exp(ltime)
ctime <- rexp(n, rate = 0.01)  #censoring_time~Exp(0.01)
simcure <- data.frame(x1 = x1,  c_time = ctime, r_time = time)
simcure$Ncure <- status
cuttime <- quantile(time, probs=0.8) 
simcure[simcure$c_time > cuttime,]$c_time <- cuttime  #admin censoring
simcure$r_time <- ifelse(simcure$Ncure == 1, simcure$r_time, Inf)
simcure$obstime <- ifelse(simcure$c_time < simcure$r_time, simcure$c_time, simcure$r_time)
simcure$event <- ifelse(simcure$obstime == simcure$c_time, 0, 1)

x.train <- simcure$x1
times <- simcure$obstime
delta <- simcure$event

#non-cured subset
notcure <- simcure$Ncure==1
tdata <- simcure[notcure,]
library(BART3)
post <- abart(x.train=tdata$x1, times=tdata$obstime, delta=tdata$event)
#plot(xb[notcure], post$yhat.train.mean)
plot(xb[notcure], ltime[notcure])
par(mfrow=c(1,3))
plot(ltime[notcure],post$yhat.train.mean,pch=20,main="abart", col = ifelse(delta[notcure]==1, 1, 2))
abline(a=0,b=1,col=2)

post1 <- qbart(x.train1=x.train, x.train2=x.train, times=times, delta=delta, q=simcure$Ncure, nskip=100, keepevery=1, ndpost=1000)
#str(post1)
plot(ltime[notcure],post1$y2hat.train.mean[notcure], pch=20, col = ifelse(delta[notcure]==1, 1, 2), main="qbart")
abline(a=0,b=1,col=2)
mean(post1$prob.train)
plot(post$yhat.train.mean, post1$y2hat.train.mean[notcure], pch=20, col = ifelse(delta[notcure]==1, 1, 2), xlab="abart", ylab="qbart", main="Non-cured log(time)")
abline(a=0,b=1,col=2)

par(mfrow=c(1,2))
plot(post1$sigma, type='l', ylab="sigma", main="trajectory of sigma")
abline(h=1,col=2)

i <- 21
plot(qnorm(post1$prob.train[,i]), type='l', ylab="p", main=paste0("trajectory of p for subj.",i))  #trajectory of p
abline(h=p,col=2)

fl <- flexsurvcure(Surv(obstime,event)~x1+meanlog(x1), data=simcure, dist="lnorm")
print(fl)

#survival curves
K <- post1$K
total.surv.train <- matrix(0,nrow=n,ncol=K)
for (i in 1:n){
    total.surv.train[i,] <- post1$surv.train.mean[(i-1)*K+(1:K)]
}
total.surv.mean <- apply(total.surv.train, 2, mean)
par(mfrow=c(1,2))
plot(c(0,post1$times),c(1,total.surv.mean),ylim=c(0,1),type='s')
points(post1$times[delta == 1],total.surv.mean[delta == 1],pch = 3,col=2)
kmfit1 <- survfit(Surv(times, delta)~1)
plot(kmfit1, mark.time=TRUE)


library(qBART)

set.seed(52120)
#simulate data; p dependent
n <- 2000  #total subj
x1 <- rnorm(n, mean = 0, sd = 4)
x2 <- sample(0:1, n, replace = TRUE)
bp <- x1*0.25 + x2 -1
p <- pnorm(bp)  #not cured
status <- rbinom(rep(1,n), rep(1,n), prob = 1-p)
b0 <- 1; b1 <- 0.5; #b2 <- 0.13
xb <- b0 + b1*x1 #+ b2*x2
ltime <- rnorm(n, mean = xb)  #log(real_time)~Normal
time <- exp(ltime)
ctime <- rexp(n, rate = 0.01)  #censoring_time~Exp(0.01)
simcure <- data.frame(x1 = x1,  c_time = ctime, r_time = time)
simcure$Ncure <- status
cuttime <- quantile(time, probs=0.8) 
simcure[simcure$c_time > cuttime,]$c_time <- cuttime  #admin censoring
simcure$r_time <- ifelse(simcure$Ncure == 1, simcure$r_time, Inf)
simcure$obstime <- ifelse(simcure$c_time < simcure$r_time, simcure$c_time, simcure$r_time)
simcure$event <- ifelse(simcure$obstime == simcure$c_time, 0, 1)

x.train <- simcure$x1
times <- simcure$obstime
delta <- simcure$event

## multiple chains
post1 <- mc.qbart(x.train1=x.train, x.train2=x.train, times=times, delta=delta, mc.cores=8)
plot(post1$sigma[,1], type='l')
plot(ltime, post1$y2hat.train.mean, pch=20)
abline(a=0,b=1,col=2)

set.seed(52520)
## Friedman function


#real data with complete covariates
library(flexsurvcure)
str(bc)
bcx <- bc$group
bct <- bc$rectime
bcd <- bc$censrec
postbc <- mc.qbart(x.train1=bcx, x.train2=bcx, times=bct, delta=bcd, mc.cores=8)
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

#real data with missing covariates
str(lung)
lx <- lung[,-(1:3)]
lt <- lung$time
ld <- lung$status-1
postl <- qbart(x.train1=lx, x.train2=lx, times=lt, delta=ld)
str(postl)
mean(postl$prob.train)

mix <- flexsurvcure(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + meal.cal + wt.loss, data=lung, dist="lnorm")
print(mix)

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

