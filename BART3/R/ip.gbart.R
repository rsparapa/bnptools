
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani
## ip.gbart

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

ip.gbart <- function(
                     x.train, y.train,
                     x.test=matrix(0,0,0), type='wbart',
                     ntype=as.integer(
                         factor(type,
                                levels=c('wbart', 'pbart', 'lbart'))),
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     lambda=NA, tau.num=c(NA, 3, 6)[ntype],
                     offset=NULL, w=rep(1, length(y.train)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     probs=c(0.025, 0.975),
                     mc.cores = 2L, nice = 19L, seed = 99L,
                     shards = 1L, weight=rep(NA, shards)
                     )
{
    if(length(x.test)==0)
        stop('Supply a non-zero length x.test matrix')

    if(shards<=2) stop('The number of shards must be >2')

    if(is.na(ntype))
        stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

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
        ##rm.const <- temp$rm.const
        rm(temp)
    }

    N = length(y.train)
    s = 2*N*(shards:1)/(shards*(shards+1))
    s[1]=s[shards]+s[1]
    s[shards]=0
    s=round(s)
    s[1]=s[1]+N-sum(s)
    strata = integer(N)
    h=0
    for(i in 1:(shards-1)) {
        strata[h+(1:s[i])]=i
        h=h+s[i]
    }

    for(h in 1:shards) {
        i=shards-h
        if(h==1) {
            Y.train = y.train[ strata==1 ]
            X.train = x.train[ , strata==1 ]
            X.test = x.train[ , strata==i ]
        } else if(h==shards) {
            Y.train = post$yhat.test.mean
            X.train = x.train[ , strata==1 ]
            X.test = x.test
        } else {
            Y.train = c(y.train[ strata==h ], post$yhat.test.mean)
            X.train = cbind(x.train[ , strata==h ], X.test)
            X.test = x.train[ , strata==i ]
        }

        post = mc.gbart(x.train=X.train,
                                  y.train=Y.train,
                                  x.test=X.test, type=type, ntype=ntype,
                                  sparse=sparse, theta=theta, omega=omega,
                                  a=a, b=b, augment=augment, rho=rho,
                                  xinfo=xinfo, usequants=usequants,
                                  rm.const=rm.const,
                                  sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                                  k=k, power=power, base=base,
                                  lambda=lambda, tau.num=tau.num,
                                  offset=offset,
                                  w=w, ntree=ntree, numcut=numcut,
                                  ndpost=ndpost, nskip=nskip,
                                  keepevery=keepevery, printevery=printevery,
                                  mc.cores=mc.cores, nice=nice, seed=seed,
                                  ##shards=shards, NO MODIFIED LISA TRICK
                                  transposed=TRUE)
        }

    return(post)
}

    ## p = nrow(x.train) ## transposed: columns of x.train
    ## q = ncol(x.test)  ## transposed: rows of x.test

    ## rs <- stratrs(y.train, shards, seed)
    ## shard. <- list()
    ## for(h in 1:shards) shard.[[h]] <- rs==h

    ## for(h in 1:shards) {
    ##     Y.train = y.train[ shard.[[h]] ]
    ##     if(h==shards)  ## (h-1)/(shards-1)=1
    ##         Y.train = post$yhat.test.mean
    ##     else if(h>1) { ## (h-1)/(shards-1)<1
    ##         rs. <- stratrs(Y.train, shards-1, seed)
    ##         Y.train[rs.<h] = post$yhat.test.mean[rs.<h]
    ##     }

    ##     if(h==shards) X.test=x.test
    ##     else X.test=x.train[ , shard.[[h+1]] ]

    ##     post = mc.gbart(x.train=x.train[ , shard.[[h]] ],

