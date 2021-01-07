
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani
## ss.gbart

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
                     x.train=matrix(0,0,0),
                     y.train=NULL,
                     x.test=matrix(0,0,0), type='wbart',
                     ntype=as.integer(
                         factor(type,
                                levels=c('wbart', 'pbart', 'lbart'))),
                     RDSfile=NULL, strata=NULL, cum.weight=TRUE,
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     varprob=NULL,
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
    if(length(y.train)==0)
        stop('Supply a non-zero length y.train vector')
    if(length(x.train)==0)
        stop('Supply a non-zero length x.train matrix')
    ## if(length(x.test)==0)
    ##     stop('Supply a non-zero length x.test matrix')

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

    N = length(y.train)
    if(length(strata)==0) strata = stratrs(y.train, shards)
    post = list()
    W=0
    for(h in 1:shards) {
        strata.h = which(strata == h)
        X.train = x.train[ , strata.h ]
        if(h==shards) X.test = x.test
        else {
            j=h+1
            strata.j = which(strata == j)
            X.test = x.train[ , strata.j ]
        }
        ## else if(h==1) X.test=X.train
        ## else {
        ##     X.test. = X.test
        ##     M = ncol(X.test.)
        ##     for(j in 1:(h-1)) {
        ##         if(j==1) X.test = X.test.[ , seq(1, M, h)]
        ##         else X.test = cbind(X.test, X.test.[ , seq(j, M, h)])
        ##     }
        ##     X.test = cbind(X.test, X.train[ , seq(h, ncol(X.train), h)])
        ## }
        if(h==1) {
            Y.train = y.train[ strata.h ]
            n = length(Y.train)
            m = 0
            w.train = rep(1, n)
            ##X.train = x.train[ , strata.h ]
            ##X.test = X.train
            treeinit=FALSE
            trees=NULL
        } else {
            Y.train = c(y.train[ strata.h ], post[[i]]$yhat.test.mean)
            if(type=='wbart') {
                ##sigdf=10
                ##sigquant=0.75
                sigest=post[[i]]$sigma.mean
            } else {
                Y.train = 1*(Y.train>0)
                ##treeinit=TRUE
                ##trees=write.trees(post[[i]]$treedraws, thin=ndpost)
            }
            n = length(Y.train)
            m = length(post[[i]]$yhat.test.mean)
            if(cum.weight)
                ##w.train = c(rep(1, n-m), rep(sqrt(m/W), m))
                w.train = c(rep(sqrt(W/m), n-m), rep(1, m))
            else {
                SD=matrix(post[[i]]$sigma., nrow=ndpost, ncol=m)
                w.train = c(rep(1, n-m),
                           sqrt(apply(post[[i]]$yhat.test/SD, 2, var)))
            }
            X.train = cbind(X.train, X.train)
            varprob=post[[i]]$varprob.mean
            ##X.train = cbind(x.train[ , strata.h ], X.test)
            ##if(h==shards) X.test = x.test
            ##else X.test = x.train[ , strata.h ]
            ##return(list(post=post[[i]], trees=trees))
        }
        W=W+n-m

        if(debug) i=h
        else i=1

        if(!file.exists(paste0(RDSfile, '.', h))) {
            print(c(shard=h))
            post[[i]] = mc.gbart(x.train=X.train, y.train=Y.train,
                                 x.test=X.test,
                                 type=type, ntype=ntype,
                                 treeinit=treeinit, trees=trees,
                                 sparse=sparse, theta=theta, omega=omega,
                                 a=a, b=b, augment=augment, rho=rho,
                                 varprob=varprob,
                                 xinfo=xinfo, usequants=usequants,
                                 rm.const=rm.const,
                                 sigest=sigest, sigdf=sigdf,
                                 sigquant=sigquant,
                                 k=k, power=power, base=base,
                                 lambda=lambda, tau.num=tau.num,
                                 offset=offset,
                                 w=w.train, ntree=ntree, numcut=numcut,
                                 ndpost=ndpost, nskip=nskip,
                                 keepevery=keepevery, printevery=printevery,
                                 mc.cores=mc.cores, nice=nice, seed=NULL,
                                 ##shards=shards, NO MODIFIED LISA TRICK
                                 transposed=TRUE)

            if(class(post[[i]])[1]!=type) return(post)
            else if(type!='wbart') post[[i]]$sigma.mean = 1

            post[[i]]$shard = h
            post[[i]]$w.train = w.train
            if(length(RDSfile)>0) saveRDS(post[[i]], paste0(RDSfile, '.', h))
        } else {
            post[[i]] = readRDS(paste0(RDSfile, '.', h))
        }
    }

    if(debug) {
        post$strata = strata
        return(post)
    } else {
        post[[1]]$strata = strata
        return(post[[1]])
    }
}

