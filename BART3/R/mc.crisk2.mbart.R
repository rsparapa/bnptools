
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


mc.crisk2.mbart <- function(
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
    id = NULL,
    seed=99,    ## mc.crisk2.mbart only
    seed2=9,    ## mc.crisk2.mbart only
    mc.cores=getOption('mc.cores', 2L), ## mc.crisk2.mbart only
    nice=19L    ## mc.crisk2.mbart only
)
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    Q = nrow(x.test)

    if(Q==0) 
        stop("'x.test' must be specified: perhaps 'x.test=x.train'")
    
    post <- mc.surv.bart(x.train=x.train,
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
                  printevery=printevery,
                  seed=seed, mc.cores=mc.cores, nice=nice)

    if(attr(post, 'class')!='survbart') return(post)

    tx.train = cbind(times, x.train)[delta>0, ]
    dimnames(tx.train) = dimnames(post$tx.train)
    
    post2 <- mc.mbart(x.train=tx.train, y.train=delta[delta>0],
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
                   printevery=printevery,
                   seed=seed2, mc.cores=mc.cores, nice=nice)

    if(attr(post2, 'class')!='mbart') return(post2)

    K <- post$K
    H <- post2$K*nrow(post$tx.test)/K ## the number of different settings
    post2$H <- H
    post2$K <- post$K
    post2$times <- post$times
    post2$surv.test <- post$surv.test
    post2$cif.test <- post2$prob.test
    post2$tx.test <- post$tx.test

    for(j in 1:K) 
        for(h in 1:H) {
            l <- H*(j-1)+h

            if(j==1) 
                post2$cif.test[ , l] <- 
                    post$prob.test[ , 1]*post2$cif.test[ , l]
            else
                post2$cif.test[ , l] <- post2$cif.test[ , l-1]+
                    post2$surv.test[ , j-1]*post$prob.test[ , j]*
                    post2$cif.test[ , l]
        }
    
    post2$cif.test.mean <- apply(post2$cif.test, 2, mean)
    post2$surv.test.mean <- apply(post2$surv.test, 2, mean)

    attr(post2, 'class') <- 'crisk2mbart'

    return(post2)
}
