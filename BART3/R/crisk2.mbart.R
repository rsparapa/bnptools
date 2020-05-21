
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


crisk2.mbart <- function(
    x.train=matrix(0,0,0), ##y.train=NULL,
    ##x.train2=x.train, y.train2=NULL,
    times=NULL, delta=NULL, K=NULL, 
    x.test=matrix(0,0,0), ##x.test2=x.test,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL, ##rho2=NULL,
    xinfo=matrix(0,0,0), ##xinfo2=matrix(0,0,0),
    usequants=FALSE,
    rm.const=TRUE, type='pbart',
    ntype=as.integer(
        factor(type, levels=c('wbart', 'pbart', 'lbart'))),
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, ##offset2 = NULL,
    tau.num=c(NA, 3, 6)[ntype],
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    printevery=100L,
    transposed=FALSE,
    id = NULL,
    seed=99,    ## mc.crisk2.mbart only
    mc.cores=2, ## mc.crisk2.mbart only
    nice=19L    ## mc.crisk2.mbart only
)
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    Q = nrow(x.test)

    if(Q==0) x.test=x.train
    
    post <- surv.bart(x.train=x.train,
                      times=times, delta=1*(delta>0), K=K,
                  x.test=x.test, type=type,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  k=k, power=power, base=base,
                  offset=offset, tau.num=tau.num,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery, 
                  printevery=printevery)

    return(post)
    if(attr(post, 'class')!='survbart') return(post)

    train = (delta>0 & post$times==post$tx.train[1, ])
    post2 <- mbart(x.train=as.matrix(post$tx.train[ , train]),
                   y.train=delta[delta>0],
                   x.test=post$tx.test, type=type,
                   sparse=sparse, theta=theta, omega=omega,
                   a=a, b=b, augment=augment, rho=rho,
                   xinfo=xinfo, usequants=usequants,
                   rm.const=rm.const,
                   k=k, power=power, base=base,
                   offset=offset, tau.num=tau.num,
                   ntree=ntree, numcut=numcut,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery, 
                   printevery=printevery, transposed=TRUE)

    if(attr(post2, 'class')!='mbart') return(post2)

    post$offset2 <- post2$offset
    post$id <- id
    post$times <- times
    post$K <- K
    post$J <- post2$K ## notation is the opposite of slides!

    if(!transposed) {
        post$tx.train <- x.train
        ##post$tx.train2 <- post2$x.train
    } else {
        post$tx.train <- t(x.train)
        ##post$tx.train2 <- t(x.train2)
    }

    post$type <- type
    post$treedraws2 <- post2$treedraws
    post$varcount2 <- post2$varcount
    post$varcount2.mean <- post2$varcount.mean
    post$varprob2 <- post2$varprob
    post$varprob2.mean <- post2$varprob.mean
    post$rm.const <- rm.const
    ##post$rm.const2 <- rm.const2
    post$yhat.train <- NULL
    post$yhat.train.mean <- NULL

    if(length(x.test)>0) {
        if(!transposed) {
            post$tx.test <- x.test
            ##post$tx.test2 <- x.test2
        } else {
            post$tx.test <- t(x.test)
            ##post$tx.test2 <- t(x.test2)
        }

        H <- nrow(post$tx.test)/K ## the number of different settings

        post$yhat.test2 <- post2$yhat.test
        post$prob.test2 <- post2$prob.test

        post$surv.test <- 1-post$prob.test
        post$cif.test <- post$prob.test*post$prob.test2
        post$cif.test2 <- post$prob.test*(1-post$prob.test2)

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$cif.test[ , l] <- post$cif.test[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test[ , l]
                post$cif.test2[ , l] <- post$cif.test2[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test2[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*
                    post$surv.test[ , l]
            }

        post$cif.test.mean <- apply(post$cif.test, 2, mean)
        post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'crisk2mbart'

    return(post)
}
