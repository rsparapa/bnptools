
library(BART3)

options(mc.cores=8)
B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

## load the advanced lung cancer example
data(lung)

N <- length(lung$status)

table(lung$ph.karno, lung$pat.karno, useNA='ifany')

## if physician's KPS unavailable, then use the patient's
h <- which(is.na(lung$ph.karno))
lung$ph.karno[h] <- lung$pat.karno[h]

times <- lung$time
delta <- lung$status-1 ##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead

## this study reports time in days rather than weeks or months
## coarsening from days to weeks or months will reduce the computational burden
##times <- ceiling(times/30)
times <- ceiling(times/7)  ## weeks

##table(times)
table(delta)

## matrix of observed covariates
##x.train <- cbind(lung$sex, lung$age, lung$ph.karno)
x.train <- data.frame(lung$sex, lung$age, lung$ph.karno, factor(lung$inst))

## lung$sex:        Male=1 Female=2
## lung$age:        Age in years
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

names(x.train) <- c('M(1):F(2)', 'age(39:82)', 'ph.karno(50:100:10)',
                    'Institution')
##dimnames(x.train)[[2]] <- c('M(1):F(2)', 'age(39:82)', 'ph.karno(50:100:10)')

table(x.train[ , 1])
summary(x.train[ , 2])
table(x.train[ , 3])

## run one long MCMC chain in one process
## set.seed(99)
## post <- surv.bart(x.train=x.train, times=times, delta=delta,
##                   x.test=x.train, sparse=TRUE, K=50, usequants=TRUE)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
                     x.test=x.train, mc.cores=B, seed=99, sparse=TRUE)

K=post$K
(NK=N*K)
x.test = rbind(post$tx.test, post$tx.test) 
## sex pushed to col 2, since time is always in col 1
x.test[ , 2]=rep(1:2, each=NK)

pred <- predict(post, x.test, mc.cores=B)

RR = matrix(nrow=post$ndpost, ncol=K)
for(j in 1:K) {
    h=seq(j, NK, K)
    RR[ , j]=apply(pred$prob.test[ , h]/pred$prob.test[ , NK+h], 1, mean)
}

RR.mean = apply(RR, 2, mean)
RR.lower = apply(RR, 2, quantile, 0.025)
RR.upper = apply(RR, 2, quantile, 0.975)
RR.mean. = apply(RR, 1, mean)
RR.lower. = quantile(RR.mean., probs=0.025)
RR.upper. = quantile(RR.mean., probs=0.975)
RR.mean. = mean(RR.mean.)

plot(post$times, RR.mean, type='l', log='y',
     ylim=c(0.1, 10), ylab='RR(t, x)', xlab='t (weeks)',
     sub="Friedman's partial dependence function")
lines(post$times, RR.lower, lty=2)
lines(post$times, RR.upper, lty=2)
abline(h=1, col=8)
abline(h=RR.mean., col=2)
abline(h=c(RR.lower., RR.upper.), col=2, lty=2)
legend('bottomleft', c('M vs. F', 'Time-varying', 'Proportionality'),
       col=c(0, 1, 2), lty=1)


