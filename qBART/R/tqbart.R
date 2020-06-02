
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

tqbart=function(x.train1=NULL, x.train2, times, delta, q,
               x.test1=matrix(0,0,0), x.test2=matrix(0,0,0), K=100,
               ## type='abart',
               ntype=1,
               sparse=FALSE, theta=0, omega=1,
               a=0.5, b=1, augment=FALSE, rho1=NULL, rho2=NULL,
               x1info=matrix(0,0,0), x2info=matrix(0,0,0),  usequants=FALSE,
               rm.const=TRUE,
               sigest=NA, sigdf=3, sigquant=0.90,
               k=2, power=2, base=0.95,
               ##sigmaf=NA,
               lambda=NA, tau.num=c(NA, 3, 6)[ntype],
               ##tau.interval=0.9973,
               ## offset=NULL,
               w=rep(1, length(times)),
               ntree=c(200L, 50L, 50L)[ntype], numcut1=100L, numcut2=100L,
               ndpost=1000L, nskip=100L,
               keepevery=c(1L, 10L, 10L)[ntype],
               printevery=100L, transposed=FALSE,
               mc.cores = 1L, nice = 19L, seed = 99L
               )
{
    if(ntype!=1) stop('ntype must be 1')

    y.train=log(times)

    n = length(y.train)

    if(n!=length(delta))
       stop("length of times and delta must be equal")

    delta=as.integer(delta)

    if(length(x.train1)==0) x.train1 = x.train2

    ## ## sort data by delta; censoring obs first
    ## if(length(x.train1)==0) x.train1 = x.train2
    ## newd <- cbind(y.train, delta, x.train1, x.train2)
    ## fulld <- newd[order(delta),]
    ## y.train = fulld[, 1]
    ## delta = fulld[, 2]
    ## x.train1 = fulld[, 3:(2+ncol(x.train1))]
    ## x.train2 = fulld[, (3+ncol(x.train1)):ncol(fulld)]
    ## cn = sum(delta == 0)  #total number of censored

    if(!transposed) {
        temp1 = bartModelMatrix(x.train1, numcut1, usequants=usequants,
                               xinfo=x1info, rm.const=rm.const)
        x.train1 = t(temp1$X)
        numcut1 = temp1$numcut
        x1info = temp1$xinfo
        temp2 = bartModelMatrix(x.train2, numcut2, usequants=usequants,
                               xinfo=x2info, rm.const=rm.const)
        x.train2 = t(temp2$X)
        numcut2 = temp2$numcut
        x2info = temp2$xinfo
        if(length(x.test1)>0) {
            x.test1 = bartModelMatrix(x.test1)
            x.test1 = t(x.test1[ , temp1$rm.const])
        }
        if(length(x.test2)>0){
            x.test2 = bartModelMatrix(x.test2)
            x.test2 = t(x.test2[ , temp2$rm.const])
        }
        rm.const1 <- temp1$rm.const
        rm.const2 <- temp2$rm.const
        rm(temp1,temp2)
    }
    else {
        rm.const1 <- rm.const2 <- NULL
    }
    
    if(n!=ncol(x.train1))
        stop('The length of times and the number of rows in x.train must be identical')
    
    if(length(x.test1)==0) x.test1 = x.test2
    p1 = nrow(x.train1); p2 = nrow(x.train2)
    np = ncol(x.test1)
    if(length(rho1)==0) rho1=p1
    if(length(rho2)==0) rho2=p2
    if(length(rm.const1)==0) rm.const1 <- 1:p1
    if(length(rm.const2)==0) rm.const2 <- 1:p2


    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side

    check1=check2=(np>0 && np==n)

    for(i in 1:n)
        for(j in 1:p1) {
            if(check1) check1=((is.na(x.train1[j, i]) && is.na(x.test1[j, i])) ||
                             (!is.na(x.train1[j, i]) && !is.na(x.test1[j, i]) &&
                              x.train1[j, i]==x.test1[j, i]))
            
            while(is.na(x.train1[j, i])) {
                h=sample.int(n, 1)
                x.train1[j, i]=x.train1[j, h]
            }
        }

    if(check1) x.test1=x.train1
    else if(np>0) {
        for(i in 1:np)
            for(j in 1:p1)
                while(is.na(x.test1[j, i])) {
                    h=sample.int(np, 1)
                    x.test1[j, i]=x.test1[j, h]
                }
    }

    for(i in 1:n)
        for(j in 1:p2) {
            if(check2) check2=((is.na(x.train2[j, i]) && is.na(x.test2[j, i])) ||
                             (!is.na(x.train2[j, i]) && !is.na(x.test2[j, i]) &&
                              x.train2[j, i]==x.test2[j, i]))
            
            while(is.na(x.train2[j, i])) {
                h=sample.int(n, 1)
                x.train2[j, i]=x.train2[j, h]
            }
        }

    if(check2) x.test2=x.train2
    else if(np>0) {
        for(i in 1:np)
            for(j in 1:p2)
                while(is.na(x.test2[j, i])) {
                    h=sample.int(np, 1)
                    x.test2[j, i]=x.test2[j, h]
                }
    }

    ## initialize offsets
    binoffset <- qnorm(mean(q))  
    offset <- mean(y.train[q==1])  #cov-adjusted center of log(times)

    y.train = y.train-offset
    
    
    if(is.na(lambda)) {
        if(is.na(sigest)) {
            if(p2 < n) {
                nqy <- q==1
                y <- y.train[nqy]
                x <- x.train1[,nqy]
                sigest = summary(lm(y~.,data.frame(t(x),y)))$sigma}
            else sigest = sd(y.train)
        }
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
    
    ptm <- proc.time()

    res = .Call("tqbart",
                ## ntype,
                ##as.integer(factor(type, levels=check))-1,
                n,  #number of observations in training data
                p1,  #dimension of x1
                p2,  #dimension of x2
                np, #number of observations in test data
                x.train1,   #pxn training data for cure status
                x.train2,   #pxn training data for y
                y.train,   #training data log(time)
                delta,     ## censoring indicator
                q,        ##initial guess of cure status
                x.test1,    #p*np test data for cure status
                x.test2,    #p*np test data for y
                ntree,
                numcut1,
                numcut2,
                ndpost*keepevery,
                nskip,
                keepevery,
                power,
                base,
                binoffset,
                offset,
                3/(k*sqrt(ntree)),
                sigdf,
                lambda,
                sigest,
                w,
                sparse,
                theta,
                omega,
                ## grp,
                a,
                b,
                rho1,
                rho2,
                augment,
                printevery,
                x1info,
                x2info
                )

    res$proc.time <- proc.time()-ptm

    res$prob.train = pnorm(res$y1hat.train)
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
    
    K <- min(n, K)
    events=unique(sort(times))
    if(length(events)>K) {
        events <- unique(quantile(times, probs=(1:K)/K))
        attr(events, 'names') <- NULL
    }
    K <- length(events)

    res$surv.train <- matrix(nrow=ndpost, ncol=n*K)
    
    for(i in 1:n)
        for(j in 1:K) {
            h <- (i-1)*K+j
            surv.ncure <- pnorm(log(events[j]),
                      mean=res$y2hat.train[ , i],
                      sd=res$sigma,
                      lower.tail=FALSE)
            res$surv.train[ , h] <- 1-res$prob.train[, i] + res$prob.train[ ,i]*surv.ncure
        }
    
    res$y2hat.train.mean <- apply(res$y2hat.train, 2, mean)
    res$surv.train.mean <- apply(res$surv.train, 2, mean)
    
    if(np>0) {
        res$prob.test = pnorm(res$y1hat.test)
        res$prob.test.mean <- apply(res$prob.test, 2, mean)
        res$surv.test <- matrix(nrow=ndpost, ncol=np*K)

        for(i in 1:np)
            for(j in 1:K) {
                h <- (i-1)*K+j
                surv.ncure <- pnorm(log(events[j]),
                          mean=res$y2hat.test[ , i],
                          sd=res$sigma[-(1:nskip)],
                          lower.tail=FALSE)
                res$surv.test[ , h] <- 1-res$prob.test[, i] + res$prob.test[, i]*surv.ncure
            }
        
        res$y2hat.test.mean <- apply(res$y2hat.test, 2, mean)
        res$surv.test.mean <- apply(res$surv.test, 2, mean)
        
    }

    res$times = events
    res$K = K
    res$binaryoffset = binoffset
    res$offset = offset
    ## names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    
    dimnames(res$varcount1)[[2]] = as.list(dimnames(x.train1)[[1]])
    dimnames(res$varprob1)[[2]] = as.list(dimnames(x.train1)[[1]])
    res$varcount1.mean <- apply(res$varcount1, 2, mean)
    res$varprob1.mean <- apply(res$varprob1, 2, mean)
    dimnames(res$varcount2)[[2]] = as.list(dimnames(x.train2)[[1]])
    dimnames(res$varprob2)[[2]] = as.list(dimnames(x.train2)[[1]])
    res$varcount2.mean <- apply(res$varcount2, 2, mean)
    res$varprob2.mean <- apply(res$varprob2, 2, mean)
    res$rm.const <- rm.const
    res$x.train1 <- x.train1
    res$x.train2 <- x.train2
    attr(res, 'class') <- 'qbart'
    return(res)
}
