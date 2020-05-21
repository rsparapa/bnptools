
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2018 Robert McCulloch and Rodney Sparapani

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

hbart=function(
               x.train, y.train,
               x.test=matrix(0,0,0), 
               sparse=FALSE, theta=0, omega=1,
               a=0.5, b=1, augment=FALSE, rho=NULL,
               xinfo=matrix(0,0,0), usequants=FALSE,
               rm.const=TRUE,
               sigest=c(NA, NA), sigdf=c(3, 3),
               sigquant=c(0.90, 0.90),
               k=c(5, 2), power=c(2, 2), base=c(0.95, 0.95),
               lambda=c(NA, NA), tau.num=c(NA, NA),
               offset=c(NA, NA), 
               ntree=c(200L, 50L), numcut=100L,
               ndpost=1000L, nskip=100L, keepevery=1L,
               printevery=100L, transposed=FALSE,
               mc.cores = 1L, nice = 19L, seed = 99L
               )
{
    n = length(y.train)

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
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

    if(n!=ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p = nrow(x.train)
    np = ncol(x.test)

    if(np>0 && p!=nrow(x.test))
        stop('The number of columns in x.train and x.test must be identical')

    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    qchi = 0
    tau = 0
    
    for(i in 1:2) { ## 1 BART, 2 HBART
        if(i==1) y. = y.train
        else y. = log(residuals(fit)^2)
        
        if(is.na(offset[i])) offset[i]=mean(y.)

        y.= y.-offset[i]

        if(p > n) {
            ## we use stepwise selection to reduce p<n
            step=srstepwise(t(x.train), y., dist='normal')
            
            if(length(step)==0) fit = lm(y.~1)
            else fit = lm(y.~I(t(x.train)[ , step]))
        }
        else fit = lm(y.~I(t(x.train)))
        
        if(is.na(lambda[i])) {
            if(is.na(sigest[i])) sigest[i] = sigma(fit)
            qchi[i] = qchisq(1-sigquant[i], sigdf[i])
            lambda[i] = (sigest[i]^2)*qchi[i]/sigdf[i] 
        } else {
            sigest[i]=sqrt(lambda[i])
        }

        if(is.na(tau.num[i])) {
            tau[i]=(max(y.)-min(y.))/(2*k[i]*sqrt(ntree[i]))
        } else {
            tau[i]=tau.num[i]/(k[i]*sqrt(ntree[i]))
        }
    }

    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side

    check=(np>0 && np==n)

    for(i in 1:n)
        for(j in 1:p) {
            if(check) check=((is.na(x.train[j, i]) && is.na(x.test[j, i])) ||
                             (!is.na(x.train[j, i]) && !is.na(x.test[j, i]) &&
                              x.train[j, i]==x.test[j, i]))

            while(is.na(x.train[j, i])) {
                h=sample.int(n, 1)
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

    res = .Call("chbart",
                n,  #number of observations in training data
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

    Y=t(matrix(y.train, nrow=n, ncol=ndpost))

    res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
    res$shat.train.mean <- apply(res$shat.train, 2, mean)
    ## SD=matrix(res$sigma[-(1:nskip)], nrow=ndpost, ncol=n)
    ## log.pdf=dnorm(Y, res$yhat.train, SD, TRUE)
    ## res$sigma.mean=mean(SD[ , 1])

    ## min.log.pdf=t(matrix(apply(log.pdf, 2, min), nrow=n, ncol=ndpost))
    ## log.CPO=log(ndpost)+min.log.pdf[1, ]-
    ##     log(apply(exp(min.log.pdf-log.pdf), 2, sum))
    ## res$LPML=sum(log.CPO)

    keeptestfits <- (np>0)

    if(keeptestfits) {
        res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
        res$shat.test.mean <- apply(res$shat.test, 2, mean)
    }

    res$offset = offset
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varcount2)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob2)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$varcount2.mean <- apply(res$varcount2, 2, mean)
    res$varprob2.mean <- apply(res$varprob2, 2, mean)
    res$rm.const <- rm.const
    res$ndpost = ndpost
    attr(res, 'class') <- 'hbart'
    return(res)
}
