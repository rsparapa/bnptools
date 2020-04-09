
library(BART3)

N = 1000   #train
NP = 200   #test
P = 5       #number of covariates
M = 8      #cores
ndpost = 1000  #draws of posterior
set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
mu = x.train[ , 1]^3
print(quantile(mu, probs=c(0.1, 0.3, 0.5, 0.7, 0.9)))
x.test=matrix(runif(NP*P, -2, 2), NP, P)
x1=c(-1, 0, 1)
x.test=cbind(c(rep(x1[1], NP), rep(x1[2], NP), rep(x1[3], NP)),
             rbind(x.test, x.test, x.test)[ , -1])

y=rnorm(N, mu)  #normal
offset=mean(y)  
T=exp(y)      #event time: lognormal
C=rexp(N, 0.05)  #censoring time: exponential
delta=(T<C)*1  #event indicator
table(delta)/N  
times=(T*delta+C*(1-delta))  #observed time
q = sample(0:1, N, replace=TRUE)  #cure status
#q = rep(1, N)
x.train = cbind(x.train,q)
new.train = x.train[which(x.train[,"q"]==1),]  #uncured subset
new.test = cbind(x.test, q=1)

post1 = mc.abart(new.train, times[q==1], delta[q==1], x.test,
                 mc.cores=M, seed=99, ndpost=ndpost)

library(qBART)
post2 = mc.qbart(x.train, times, delta, q, new.test,
                 mc.cores=M, seed=99, ndpost=ndpost)


Z=2

plot(post2$yhat.test.mean, post1$yhat.test.mean, asp=1)
abline(a=0, b=1)

## plot(mu, post2$yhat.train.mean, asp=1,
##      xlim=c(-Z, Z), ylim=c(-Z, Z))
## abline(a=0, b=1)

## plot(post$yhat.train.mean, post2$yhat.train.mean, asp=1,
##      xlim=c(-Z, Z), ylim=c(-Z, Z))
## abline(a=0, b=1)

K <- post$K
par(mfrow=c(3, 1))
for(i in 1:length(x1)) {
        plot(c(0, post$times),
             c(1, pnorm(log(post$times), mean=x1[i]^3,
                   lower.tail=FALSE)),
             type='l', ylim=0:1, xlab='t', ylab='S(t, x)')
        lines(c(0, post$times),
              c(1, post$surv.test.mean[(i-1)*K+1:K]),
              col=2, type='s')
        post$surv.test.025 <- matrix(nrow=ndpost, ncol=K)
        post$surv.test.975 <- matrix(nrow=ndpost, ncol=K)
        for(j in 1:K) {
            post$surv.test.025[ , j] <-
                apply(post$surv.test[ , (i-1)*K*NP+seq(j, K*NP, K)],
                      1, mean)
            post$surv.test.975[ , j] <-
                apply(post$surv.test[ , (i-1)*K*NP+seq(j, K*NP, K)],
                      1, mean)
        }
        post$surv.test.025 <- apply(post$surv.test.025, 2,
                                    quantile, probs=0.025)
        post$surv.test.975 <- apply(post$surv.test.975, 2,
                                    quantile, probs=0.975)
        lines(c(0, post$times), c(1, post$surv.test.025),
              col=2, type='s', lty=2)
        lines(c(0, post$times), c(1, post$surv.test.975),
              col=2, type='s', lty=2)
}
par(mfrow=c(1, 1))
