
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani
## s.gbart

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

ss.gbart <- function(
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
                     offset=NULL, ##w=rep(1, length(y.train)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     probs=c(0.025, 0.975),
                     mc.cores = 2L, nice = 19L, seed = 99L,
                     shards = 1L, weight=rep(NA, shards),
                     debug=FALSE
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

    if(length(offset)==0) {
        offset=mean(y.train)
        if(type=='pbart') offset=qnorm(offset)
        else if(type=='lbart') offset=qlogis(offset)
    }

    if(type=='wbart') y.train = y.train-offset

    N = length(y.train)
    strata = stratrs(y.train, shards)
    post = list() ## DEBUGGING
    ## pbart=(type=='pbart')
    W=0
    for(h in 1:shards) {
        strata.h = which(strata == h)
        if(h==1) {
            z.train = y.train[ strata.h ]
            Y.train = z.train
            n = length(Y.train)
            m = 0
            w.train = rep(1, n)
            X.train = x.train[ , strata.h ]
            X.test = X.train
        } else {
            z.train = c(y.train[ strata.h ], post[[h-1]]$yhat.test.mean)
            Y.train = z.train
            if(type!='wbart') Y.train = 1*(z.train>0)
            n = length(Y.train)
            m = length(post[[h-1]]$yhat.test.mean)
            w.train = c(rep(1, n-m), rep(sqrt(m/W), m))
            X.train = cbind(x.train[ , strata.h ], X.test)
            if(h==shards) X.test = x.test
            else X.test = x.train[ , strata.h ]
        }
        W=W+n-m
        print(c(h=h, W=W, n-m))

        post[[h]] = mc.gbart(x.train=X.train, y.train=Y.train, x.test=X.test,
                        z.train=z.train, type=type, ntype=ntype,
                                  sparse=sparse, theta=theta, omega=omega,
                                  a=a, b=b, augment=augment, rho=rho,
                                  xinfo=xinfo, usequants=usequants,
                                  rm.const=rm.const,
                                  sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                                  k=k, power=power, base=base,
                                  lambda=lambda, tau.num=tau.num,
                                  offset=offset,
                                  w=w.train, ntree=ntree, numcut=numcut,
                                  ndpost=ndpost, nskip=nskip,
                                  keepevery=keepevery, printevery=printevery,
                                  mc.cores=mc.cores, nice=nice, seed=seed,
                                  ##shards=shards, NO MODIFIED LISA TRICK
                                  transposed=TRUE)

        if(class(post[[h]])[1]!=type) return(post)
        else if(type!='wbart') post[[h]]$sigma.mean = 1
        ## else if(type=='wbart' && pbart) {
        ##     class(post[[h]])='pbart'
        ##     post[[h]]$prob.test=pnorm(post[[h]]$yhat.test)
        ##     post[[h]]$prob.test.mean=apply(post[[h]]$prob.test, 2, mean)
        ## }
    }

    if(debug) {
        post$strata = strata
        return(post)
    } else {
        post[[h]]$strata = strata
        return(post[[h]])
    }
}

