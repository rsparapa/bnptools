##################################################
## This is just a stub (runs fast) example for testing.
##  For more realistic examples, please see:
##   (i) the vignette at www.rob-mcculloch.org
##   (ii) the example simulated data (see ?simdat)
##        and the longer run in ?hbartonsimd, 
##        where a saved run of hbart is run on simdat is plotted. 
##################################################
library(hbart)
##simulate data
set.seed(99)

# train data
n=500 #train data sample size
p=10 
x = matrix(runif(n*p),ncol=p) #iid uniform x values
fx = 4*(x[,1]^2) #quadratric function f
sx = .2*exp(2*x[,1]) # exponential function s
y = fx + sx*rnorm(n) # y = f(x) + s(x) Z

#test data (the p added to the variable names is for predict)
np=500 #test data sample size
xp = matrix(runif(np*p),ncol=p)
i=order(xp[,1])
xp=xp[i,]
fxp = 4*(xp[,1]^2)
sxp = .2*exp(2*xp[,1])
yp = fxp + sxp*rnorm(np)

##run hbart MCMC
# The number of interations is kept small to make example run,
##!!!!  REAL APPLICATIONS MAY NEED LONGER RUNS !!!!
#   nskip: burn in draws,
#   ndpost:kept draws,
#   nadapt: initial draws to tune MCMC,
#   numcut: number of cutpoints used for each x
#   k: bigger k gives smoother f (default is 2)
set.seed(19)
res = hbart(x,y,nskip=10,ndpost=20,nadapt=0,numcut=1000,k=5,
            summarystats=TRUE) #again, this is way too short a run!!!
## now predict to get inference
resp = predict(res,x.test=xp)

##check out of sample fit
cat("out of sample cor(f,fhat) is ",cor(fxp,resp$mmean),"\n")
cat("out of sample cor(s,shat) is ",cor(sxp,resp$smean),"\n")

##plot estimated vs. true
##plot the data
plot(xp[,1],yp,cex.axis=1.5,cex.lab=1.5)
lines(xp[,1],fxp,col="blue")
lines(xp[,1],fxp+2*sxp,col="blue",lty=2)
lines(xp[,1],fxp-2*sxp,col="blue",lty=2)

## add the fit
lines(xp[,1],resp$mmean) #estimate of f
lines(xp[,1],resp$mmean+2*resp$smean) #estimate of sd
lines(xp[,1],resp$mmean-2*resp$smean) #estimate of sd

print(res$mu.varprob*p)
print(res$sd.varprob*p)
