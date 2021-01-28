
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani
## ml.gbart

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

ml.gbart <- function(
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
                     shards = 1L, weight=rep(NA, shards),
                     meta=FALSE
                     )
{
    if(length(x.test)==0)
        stop('Supply a non-zero length x.test matrix')

    if(shards<=1) stop('The number of shards must be >1')

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

    rs <- stratrs(y.train, shards, seed)

    post.list = list()
    shard. <- as.list(1:shards)
    p = nrow(x.train) ## transposed: columns of x.train
    q = ncol(x.test)  ## transposed: rows of x.test

    if(meta && is.na(weight[1])) weight=rep(1, shards)

    for(h in 1:shards) {
        shard.[[h]] <- rs==h
        post.list[[h]] = mc.gbart(x.train=x.train[ , shard.[[h]] ],
                                  y.train=y.train[ shard.[[h]] ],
                                  x.test=x.test, type=type, ntype=ntype,
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
                                  shards=shards, meta=meta, transposed=TRUE)

        if(h==1) {
            post <- post.list[[1]]

            if(attr(post, 'class')!=type) return(post)

            post$weight=double(shards)
            post$treedraws=list()
            post$varcount=array(dim=c(shards, ndpost, p))
            post$varprob=array(dim=c(shards, ndpost, p))
            post$varcount.mean=matrix(nrow=shards, ncol=p)
            post$varprob.mean=matrix(nrow=shards, ncol=p)
            post$yhat.train <- NULL
            if(type=='wbart') {
                post$yhat.train.mean <- NULL
                post$yhat.test.mean <- NULL
                post$yhat.test.mean <- NULL
            } else {
                post$prob.train <- NULL
                post$prob.train.mean <- NULL
                post$prob.test.mean <- NULL
            }
        } else {
            post$offset[h] = post.list[[h]]$offset
        }

        if(type=='wbart') {
            if(h==1) post$sigma = cbind(apply(post$sigma, 1, mean))
            else post$sigma = cbind(post$sigma, apply(post.list[[h]]$sigma, 1, mean))

            if(is.na(weight[h]))
                post$weight[h] = mean(apply((post.list[[h]]$sigma[-(1:nskip), ])^(-2), 1, mean))
            else post$weight[h] = weight[h]
            ## post$weight[[h]]=matrix(apply((post.list[[h]]$sigma[-(1:nskip), ])^2, 1, mean),
            ##                         nrow=ndpost, ncol=q)
            ## post$weight[[h]]=matrix(1/apply(post$weight[[h]]*post$yhat.test, 2, var),
            ##                         nrow=ndpost, ncol=q, byrow=TRUE)
        } else {
            if(is.na(weight[h]))
                post$weight[h] = 1
            else post$weight[h] = weight[h]
            ## post$weight[[h]]=matrix(1/apply(post$yhat.test, 2, var),
            ##                         nrow=ndpost, ncol=q, byrow=TRUE)
        }

        post$treedraws[[h]]=post.list[[h]]$treedraws
        post$varcount[h, , ]=post.list[[h]]$varcount
        post$varprob[h, , ]=post.list[[h]]$varprob
        ## post$varcount[[h]]=post.list[[h]]$varcount
        ## post$varprob[[h]]=post.list[[h]]$varprob
        post$varcount.mean[h, ]=post.list[[h]]$varcount.mean
        post$varprob.mean[h, ]=post.list[[h]]$varprob.mean

        ## if(h==1) {
        ##     post$varcount.mean=apply(post$varcount[[h]], 2, mean)
        ##     post$varprob.mean=apply(post$varprob[[h]], 2, mean)/shards
        ## } else {
        ##     post$varcount.mean=post$varcount.mean+
        ##         apply(post$varcount[[h]], 2, mean)
        ##     post$varprob.mean=post$varprob.mean+
        ##         apply(post$varprob[[h]], 2, mean)/shards
        ## }
    }

    post$weight=post$weight/sum(post$weight)

    for(h in 1:shards) {
        ## if(h==1) post$yhat.test = post$weight[[1]]*post$yhat.test/post$weight.
        ## else post$yhat.test = post$yhat.test+post$weight[[h]]*post.list[[h]]$yhat.test/post$weight.
        if(h==1) {
            post$yhat.test.var = apply(post$yhat.test, 2, var)
            if(type!='wbart') post$prob.test.var = apply(post$prob.test, 2, var)
            ## if(type=='wbart') post$yhat.test.var = apply(post$yhat.test, 2, var)
            ## else post$prob.test.var = apply(post$prob.test, 2, var)
            post$yhat.test = post$weight[1]*(post$yhat.test-post$offset[1])
        } else {
            post$yhat.test.var = rbind(post$yhat.test.var,
                                       apply(post.list[[h]]$yhat.test, 2, var))
            if(type!='wbart') post$prob.test.var = rbind(post$prob.test.var,
                                       apply(post.list[[h]]$prob.test, 2, var))
            ## if(type=='wbart') post$yhat.test.var = rbind(post$yhat.test.var,
            ##                            apply(post.list[[h]]$yhat.test, 2, var))
            ## else post$prob.test.var = rbind(post$prob.test.var,
            ##                            apply(post.list[[h]]$prob.test, 2, var))
            post$yhat.test = post$yhat.test+
                 post$weight[h]*(post.list[[h]]$yhat.test-post$offset[h])
        }
    }

    post$yhat.test = post$yhat.test+mean(post$offset)

##    if(type=='wbart') {
        post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
        post$yhat.test.sd <- sqrt(apply(post$yhat.test.var, 2, mean)/shards)
        if(meta) {
            post$yhat.test.lower <- post$yhat.test.mean +
                qnorm(min(probs))*post$yhat.test.sd
            post$yhat.test.upper <- post$yhat.test.mean +
                qnorm(max(probs))*post$yhat.test.sd
        } else {
            post$yhat.test.lower <- apply(post$yhat.test, 2, quantile,
                                          probs=min(probs))
            post$yhat.test.upper <- apply(post$yhat.test, 2, quantile,
                                          probs=max(probs))
        }
##} else {
    if(type!='wbart') {
        if(type=='pbart') post$prob.test <- pnorm(post$yhat.test)
        else if(type=='lbart') post$prob.test <- plogis(post$yhat.test)
        post$prob.test.mean <- apply(post$prob.test, 2, mean)
        post$prob.test.sd <- sqrt(apply(post$prob.test.var, 2, mean)/shards)
        if(meta) {
            post$prob.test.lower <- post$prob.test.mean +
                qnorm(min(probs))*post$prob.test.sd
            post$prob.test.upper <- post$prob.test.mean +
                qnorm(max(probs))*post$prob.test.sd
        } else {
            post$prob.test.lower <- apply(post$prob.test, 2, quantile,
                                          probs=min(probs))
            post$prob.test.upper <- apply(post$prob.test, 2, quantile,
                                          probs=max(probs))
        }
    }

    return(post)
}

## old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p))
## old.stop <- nchar(old.text)

## post$treedraws$trees <- sub(old.text,
##                             paste0(as.character(shards*mc.cores*mc.ndpost), ' ',
##                                    as.character(ntree), ' ',
##                                    as.character(p)),
##                             post$treedraws$trees)

## keeptestfits <- length(x.test)>0

## for(h in 1:shards)
##     for(i in 1:mc.cores)
##         if(!(h==1 && i==1)) {
##             post$yhat.shard[(i-1)*mc.ndpost+1:mc.ndpost, shard.[[h]] ] <- post.list[[h]][[i]]$yhat.train

##             post$varcount <- rbind(post$varcount, post.list[[h]][[i]]$varcount)

##             post$treedraws$trees <- paste0(post$treedraws$trees,
##                                            substr(post.list[[h]][[i]]$treedraws$trees, old.stop+2,
##                                                   nchar(post.list[[h]][[i]]$treedraws$trees)))
##         }

## attr(post, 'class') <- type

## if(keeptestfits) {
##     post$x.test <- bartModelMatrix(x.test)
##     post$yhat.test <- predict(post, newdata=post$x.test, mc.cores=mc.cores)
##     if(type=='wbart')
##         post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
##     else {
##         if(type=='pbart')
##             post$prob.test <- pnorm(post$yhat.test)
##         else if(type=='lbart')
##             post$prob.test <- plogis(post$yhat.test)
##         post$prob.test.mean <- apply(post$prob.test, 2, mean)
##     }
## } else {
##     post$x.train <- t(x.train)
##     post$yhat.train <- predict(post, newdata=post$x.train, mc.cores=mc.cores)
##     if(type=='wbart')
##         post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
##     else {
##         if(type=='pbart')
##             post$prob.train <- pnorm(post$yhat.train)
##         else if(type=='lbart')
##             post$prob.train <- plogis(post$yhat.train)
##         post$prob.train.mean <- apply(post$prob.train, 2, mean)
##     }
## }

## return(post)
## }
