
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

liobart=function(
               x.train, y.train, x.test=matrix(0,0,0),
               sparse=FALSE, theta=0, omega=1,
               a=0.5, b=1, augment=FALSE, rho=NULL,
               N = length(y.train), ## for convenience, not an argument
               states=c(N, rep(0, N-1)), C=rep(1, N),
               neal.m=2L, m0=0, k0=0.2, k0.a=1.5, k0.b=7.5, k0.draw=1,
               a0=1.5, b0=0.5, b0.a=0.5, b0.b=1, b0.draw=1,
               alpha=1, alpha.a=0.1, alpha.b=0.1, alpha.draw=1,
               standard=FALSE,
               xinfo=matrix(0,0,0), usequants=FALSE,
               rm.const=TRUE,
               sigest=NA, sigdf=3, sigquant=0.90,
               k=2, power=2, base=0.95,
               lambda=NA, tau.num=NA,
               offset=NULL,
               ntree=50L, numcut=100L,
               ndpost=2000L, nskip=1000L,
               keepevery=10L,
               printevery=100L, transposed=FALSE,
               probs=c(0.025, 0.975),
               mc.cores = 1L, nice = 19L, seed = 99L
               )
{
    stopifnot(neal.m>=0 & neal.m==as.integer(neal.m))

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if(length(x.test)>0) {
            x.test = t(bartModelMatrix(x.test))
            if(class(rm.const)[1]=='logical' && rm.const)
                x.test = rbind(x.test[temp$rm.const, ])
        }
        rm.const <- temp$rm.const
        grp <- temp$grp
        rm(temp)
    }
    else {
        rm.const <- NULL
        grp <- NULL
    }

    if(N!=ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p = nrow(x.train)
    np = ncol(x.test)

    if(np>0 && p!=nrow(x.test))
        stop('The number of columns in x.train and x.test must be identical')

    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    if(length(offset)==0) offset=mean(y.train)

    y.train = y.train-offset

    if(p > N) {
        ## we use stepwise selection to reduce p<n
        step=srstepwise(t(x.train), y.train, dist='normal')

        if(length(step)==0) fit = lm(y.train~1)
        else fit = lm(y.train~I(t(x.train)[ , step]))
    }
    else fit = lm(y.train~I(t(x.train)))

    if(is.na(lambda)) {
        if(is.na(sigest)) sigest = sigma(fit)
        qchi = qchisq(1-sigquant, sigdf)
        lambda = (sigest^2)*qchi/sigdf #lambda parameter for sigma prior
    } else {
        sigest=sqrt(lambda)
    }

    if(is.na(tau.num)) {
        tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
    } else {
        tau=tau.num/(k*sqrt(ntree))
    }

    if(!standard) {
        scale2=sigest^2
        m0=m0+mean(residuals(fit))
        a0=a0+0.5
        k0.b=k0.b/scale2
        b0.b=b0.b/scale2
    } else {
        scale2 = 1
    }

    phi=matrix(c(m0, 1/scale2), nrow=N+neal.m, ncol=2, byrow=TRUE)
    dimnames(phi)[[2]]=c('mu', 'tau')
    prior=list(m=as.integer(neal.m),
               m0=m0, k0.a=k0.a, k0.b=k0.b,
               a0=a0, b0.a=b0.a, b0.b=b0.b,
               alpha.a=alpha.a, alpha.b=alpha.b)
    hyper=list(alpha=alpha, alpha.draw=alpha.draw,
               k0=k0, b0=b0, k0.draw=k0.draw, b0.draw=b0.draw)
    C=as.integer(C-1) ## convert from R to C/C++ indexing
    states=as.integer(states)

    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side

    check=(np>0 && np==N)

    for(i in 1:N)
        for(j in 1:p) {
            if(check) check=((is.na(x.train[j, i]) && is.na(x.test[j, i])) ||
                             (!is.na(x.train[j, i]) && !is.na(x.test[j, i]) &&
                              x.train[j, i]==x.test[j, i]))

            while(is.na(x.train[j, i])) {
                h=sample.int(N, 1)
                x.train[j, i]=x.train[j, h]
            }
        }

    if(check) x.test=x.train
    else if(np>0) {
        for(i in 1:np)
            for(j in 1:p)
                while(is.na(x.test[j, i])) {
                    h=sample.int(np, 1)
                    x.test[j, i]=x.test[j, h]
                }
    }

    ptm <- proc.time()

    ## return(list(phi, C, states, prior, hyper,
    ##             N,  #number of observations in training data
    ##             p,  #dimension of x
    ##             np, #number of observations in test data
    ##             x.train,   #pxn training data x
    ##             y.train,   #pxn training data x
    ##             x.test,    #p*np test data x
    ##             ntree,
    ##             numcut,
    ##             ndpost*keepevery,
    ##             nskip,
    ##             keepevery,
    ##             power,
    ##             base,
    ##             offset,
    ##             tau,
    ##             sigdf,
    ##             lambda,
    ##             sigest,
    ##             ##w,
    ##             sparse,
    ##             theta,
    ##             omega,
    ##             grp,
    ##             a,
    ##             b,
    ##             rho,
    ##             augment,
    ##             printevery,
    ##             xinfo
    ##             ))

    res = .Call("cliobart",
                phi, C, states, prior, hyper,
                N,  #number of observations in training data
                p,  #dimension of x
                np, #number of observations in test data
                x.train,   #pxn training data x
                y.train,   #pxn training data x
                x.test,    #p*np test data x
                ntree,
                numcut,
                ndpost*keepevery,
                nskip,
                keepevery,
                power,
                base,
                offset,
                tau,
                sigdf,
                lambda,
                sigest,
                ##w,
                sparse,
                theta,
                omega,
                grp,
                a,
                b,
                rho,
                augment,
                printevery,
                xinfo
                )

    res$proc.time <- proc.time()-ptm
##    res$hostname <- hostname

    Y=t(matrix(y.train, nrow=N, ncol=ndpost))

        res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
        res$yhat.train.lower <- apply(res$yhat.train, 2, quantile,
                                      probs=min(probs))
        res$yhat.train.upper <- apply(res$yhat.train, 2, quantile,
                                      probs=max(probs))
        SD=matrix(res$sigma[-(1:nskip)], nrow=ndpost, ncol=n)
        ##CPO=1/apply(1/dnorm(Y, res$yhat.train, SD), 2, mean)
        log.pdf=dnorm(Y, res$yhat.train, SD, TRUE)
        res$sigma.mean=mean(SD[ , 1])

    min.log.pdf=t(matrix(apply(log.pdf, 2, min), nrow=n, ncol=ndpost))
    log.CPO=log(ndpost)+min.log.pdf[1, ]-
        log(apply(exp(min.log.pdf-log.pdf), 2, sum))
    res$LPML=sum(log.CPO)

    keeptestfits <- (np>0)

    if(keeptestfits) {
            res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
            res$yhat.test.lower <- apply(res$yhat.test, 2, quantile,
                                         probs=min(probs))
            res$yhat.test.upper <- apply(res$yhat.test, 2, quantile,
                                         probs=max(probs))
    }

    res$offset = offset
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    res$ndpost = ndpost
    attr(res, 'class') <- 'liobart'
    return(res)
}
