
library(hbart)
##simulate data
set.seed(99)

n=500 #train data sample size
p=100 
x = matrix(runif(n*p),ncol=p) #iid uniform x values
fx = 4*(x[,1]^2) #quadratric function f
sx = .2*exp(2*x[,1]) # exponential function s
y = fx + sx*rnorm(n) # y = f(x) + s(x) Z

np=500 #test data sample size
xp = matrix(runif(np*p),ncol=p)
i=order(xp[,1])
xp=xp[i,]
fxp = 4*(xp[,1]^2)
sxp = .2*exp(2*xp[,1])
yp = fxp + sxp*rnorm(np)

##run hbart MCMC
rerun=TRUE
file.='res.rds'
if(file.exists(file.) && !rerun) {
    res=readRDS(file.)
    XPtr=FALSE
} else {
    set.seed(19)
    res = hbart(x,y,tc=8)
    ##     ##nskip=10,ndpost=20,nadapt=10)
    saveRDS(res, file.)
    XPtr=TRUE
}
## now predict to get inference
resp = predict(res,x.test=xp,XPtr=XPtr,tc=8)

if(XPtr) {

    ##check out of sample fit
    cat("out of sample cor(f,fhat) is ",cor(fxp,resp$mmean)^2,"\n")
    print(cor(resp$mmean, resp$f.test.mean)^2)
    plot(fxp, resp$f.test.mean, pch='.', col=2, asp=1)
    points(fxp, resp$mmean, pch='.', col=4)
    abline(a=0, b=1)

    cat("out of sample cor(s,shat) is ",cor(sxp,resp$smean)^2,"\n")
    print(cor(resp$smean, resp$s.test.mean)^2)
    plot(sxp, resp$s.test.mean, pch='.', col=2, asp=1)
    points(sxp, resp$smean, pch='.', col=4)
    abline(a=0, b=1)

    plot(resp$smean, resp$s.test.mean, pch='.', asp=1)
    abline(a=0, b=1)

    if(FALSE) {
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

        plot(res$mu.varprob, ylim=0:1,
             pch=1+44*(res$mu.varprob<1/p))
        points(res$sd.varprob, col=2,
               pch=1+44*(res$sd.varprob<1/p))
        abline(h=c(0, 1/p), v=1.5)
        legend('topright', c('mu', 'sd'), col=1:2, pch=1)
    }
}

