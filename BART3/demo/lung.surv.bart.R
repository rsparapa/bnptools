
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

x.test = post$tx.test[1:2, ]
## sex pushed to col 2, since time is always in col 1
x.test[ , 2]=1:2
pred. <- FPD(post, x.test, S=2, mc.cores=B)

K=post$K
plot(c(0, post$times), c(1, pred.$surv.test.mean[1:K]), type='s', col='blue',
     ylim=0:1, ylab='S(t, x)', xlab='t (weeks)',
     sub="Friedman's partial dependence function")
lines(c(0, post$times), c(1, pred.$surv.test.lower[1:K]), col='blue', type='s', lty=2)
lines(c(0, post$times), c(1, pred.$surv.test.upper[1:K]), col='blue', type='s', lty=2)
lines(c(0, post$times), c(1, pred.$surv.test.mean[K+1:K]), col='red', type='s')
lines(c(0, post$times), c(1, pred.$surv.test.lower[K+1:K]), col='red', type='s', lty=2)
lines(c(0, post$times), c(1, pred.$surv.test.upper[K+1:K]), col='red', type='s', lty=2)

##pred2=SHAP(post, pre$tx.train, tx.test, 1:2)

## plot(c(0, pre$times), c(1, pd.mu[males]), type='s', col='blue',
##      ylim=0:1, ylab='S(t, x)', xlab='t (weeks)',
##      sub="Friedman's partial dependence function")
## lines(c(0, pre$times), c(1, pd.025[males]), col='blue', type='s', lty=2)
## lines(c(0, pre$times), c(1, pd.975[males]), col='blue', type='s', lty=2)
## lines(c(0, pre$times), c(1, pd.mu[females]), col='red', type='s')
## lines(c(0, pre$times), c(1, pd.025[females]), col='red', type='s', lty=2)
## lines(c(0, pre$times), c(1, pd.975[females]), col='red', type='s', lty=2)
## if(figures!='NONE')
##     dev.copy2pdf(file=paste(figures, 'lung.pdf', sep='/'))

