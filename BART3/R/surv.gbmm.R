
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

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


surv.gbmm <- function(
    x.train = matrix(0,0,0),
    y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0,0,0),
    K=NULL, events=NULL, ztimes=NULL, zdelta=NULL,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0,0,0), usequants=FALSE,
    ## cont=FALSE,
    rm.const=TRUE, type='pbart', ntype=2,
    u.train=NULL, B=NULL,
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, tau.num=c(NA, 3, 6)[ntype],
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L, keepevery = 10L,
    printevery=100L,
    id = NULL,     ## only used by surv.bart
    seed = 99L,    ## only used by mc.surv.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(type!='pbart' || ntype!=2)
        stop("type argument must be set to 'pbart'")

    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(rho)==0) rho=ncol(x.train)

    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test, K=K,
                             events=events,
                             ztimes=ztimes, zdelta=zdelta,
                             u.train=u.train)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test

        times   <- pre$times
        K       <- pre$K
        u.train <- pre$u.train
    }
    else {
        if(length(unique(sort(y.train)))>2)
            stop('y.train has >2 values; make sure you specify times=times & delta=delta')

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    post <- gbmm(x.train=x.train, y.train=y.train,
                  x.test=x.test, type=type,
                  u.train=u.train, B=B,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  k=k, power=power, base=base,
                  offset=offset, tau.num=tau.num,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip, keepevery=keepevery,
                  printevery=printevery)

    if(type!=attr(post, 'class')) return(post)

    post$id <- id
    post$times <- times
    post$K <- K
    post$tx.train <- x.train
    post$type <- type

    ## if(keeptrainfits) {
    ##      post$surv.train <- 1-post$prob.train

    ##      H <- nrow(x.train)/K ## the number of different settings

    ##      for(h in 1:H)
    ##          for(j in 2:K) {
    ##              l <- K*(h-1)+j

    ##              post$surv.train[ , l] <-
    ## post$surv.train[ , l-1]*post$surv.train[ , l]
    ##          }

    ##      post$surv.train.mean <- apply(post$surv.train, 2, mean)
    ## }

    if(length(x.test)>0) {
        post$tx.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        post$surv.test <- 1-post$prob.test

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$surv.test[ , l] <-
                    post$surv.test[ , l-1]*post$surv.test[ , l]
            }

        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'survgbmm'

    return(post)
}
