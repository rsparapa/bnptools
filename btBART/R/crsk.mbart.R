
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

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


crsk.mbart <- function(
    x.train=matrix(0.0, 0L, 0L), y.train=NULL,
    times=NULL, delta=NULL, K=NULL,
    x.test=matrix(0.0, 0L, 0L), 
    sparse=FALSE, a=0.5, b=1, augment=FALSE, rho=NULL, 
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    cont=FALSE, rm.const=TRUE, 
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = 0,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    nkeeptrain=ndpost, nkeeptest=ndpost,
    ##nkeeptestmean=ndpost,
    nkeeptreedraws=ndpost,
    printevery=100L,
    ##treesaslists=FALSE,
    keeptrainfits=TRUE,
    id = NULL,
    seed=99,    ## mc.crsk.mbart only
    mc.cores=2, ## mc.crsk.mbart only
    nice=19L    ## mc.crsk.mbart only
)
{

    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(rho)==0) rho=ncol(x.train)

    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test, K=K)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test

        times   <- pre$times
        K       <- pre$K
    }
    else {
        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    if(length(xinfo)==0) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               cont=cont, xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        rm.const <- temp$rm.const
        rm(temp)

        transposed <- TRUE
    }
    else {
        rm.const <- 1:ncol(x.train)
        transposed <- FALSE
    }

    post <- mbart(x.train=x.train, y.train=y.train, x.test=x.test,
                 sparse=sparse, a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  cont=cont, rm.const=rm.const,
                  k=k, power=power, base=base,
                  binaryOffset=binaryOffset,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery, nkeeptrain=0,
                  nkeeptest=nkeeptest, #nkeeptestmean=nkeeptestmean,
                  nkeeptreedraws=nkeeptreedraws, printevery=printevery,
                  transposed=transposed) #, treesaslists=treesaslists)

    if('mbart'!=attr(post, 'class')) return(post)

    C <- post$C

    if(C!=3) 
        stop('C must be 3: no event (position 1), cause 1 (2) and cause 2 (3)')

    post$binaryOffset <- binaryOffset
    post$id <- id
    post$times <- times
    post$K <- K

    if(!transposed) {
        post$tx.train <- x.train
    } else {
        post$tx.train <- t(x.train)
    }

    post$rm.const <- rm.const
    post$yhat.train <- NULL
    post$yhat.train.mean <- NULL

    if(length(x.test)>0) {
        if(!transposed) {
            post$tx.test <- x.test
        } else {
            post$tx.test <- t(x.test)
        }

        H <- nrow(post$tx.test)/K ## the number of different settings
        
        h <- seq(1, H*K*C, C) ## no event: position 1
        
        post$surv.test <- post$prob.test[ , h]
        
        post$cif.test <- post$prob.test[ , h+1]
        post$cif.test2 <- post$prob.test[ , h+2]

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$cif.test[ , l] <- post$cif.test[ , l-1]+post$surv.test[ , l-1]*post$cif.test[ , l]
                post$cif.test2[ , l] <- post$cif.test2[ , l-1]+post$surv.test[ , l-1]*post$cif.test2[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*post$surv.test[ , l]
            }

        post$cif.test.mean <- apply(post$cif.test, 2, mean)
        post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'crskmbart'

    return(post)
}
