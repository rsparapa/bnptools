
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


mc.crisk3.bart <- function(
    x.train=matrix(0,0,0), x.train2=x.train, x.train3=x.train,
    y.train=NULL, y.train2=NULL, y.train3=NULL,
    times=NULL, delta=NULL, K=NULL, events=NULL,
    x.test=matrix(0,0,0), x.test2=x.test, x.test3=x.test,
    cond=NULL, cond2=NULL,
    sparse=FALSE, theta=0, omega=1, a=0.5, b=1, augment=FALSE,
    rho=NULL, rho2=NULL, rho3=NULL,
    xinfo=matrix(0,0,0), xinfo2=matrix(0,0,0), xinfo3=matrix(0,0,0),
    usequants=FALSE, rm.const=TRUE, type='pbart',
    ntype=as.integer(
        factor(type, levels=c('wbart', 'pbart', 'lbart'))),
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, offset2 = NULL, offset3 = NULL,
    tau.num=c(NA, 3, 6)[ntype],
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L, printevery=100L,
    ## keeptestfits = NULL,
    id = NULL,
    seed=99,    ## mc.crisk3.bart only
    mc.cores=2, ## mc.crisk3.bart only
    nice=19L    ## mc.crisk3.bart only
)
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    ## handling rm.const bug by turning it off unless explicitly set
    ## if(class(rm.const)[1]=='logical') rm.const <- FALSE

    H <- 1
    Mx <- 2^31-1

    if(length(y.train)==0) {
        pre <- crisk3.pre.bart(times, delta, x.train, x.test,
                              x.train2, x.test2, x.train3, x.test3,
                              K=K, events=events)
        ## if(class(rm.const)[1]=='logical' && !rm.const) p <- ncol(pre$tx.train)
        ## else p <- ncol(pre$tx.train[ , rm.const])
        ## p <- ncol(pre$tx.train)
        Nx <- 2*max(nrow(pre$tx.train), nrow(pre$tx.test))
        keeptestfits <- (length(pre$tx.test)>0)
    } else {
        ## if(class(rm.const)[1]=='logical' && !rm.const) p <- ncol(x.train)
        ## else p <- ncol(x.train[ , rm.const])
        ## p <- ncol(x.train)
        Nx <- 2*max(nrow(x.train), nrow(x.test))
        keeptestfits <- (length(x.test)>0)
    }

    if(Nx>Mx%/%ndpost) {
        H <- ceiling(ndpost / (Mx %/% Nx))
        ndpost <- ndpost %/% H
        ##nrow*ndpost>2Gi: due to the 2Gi limit in sendMaster
        ##(unless this limit was increased): reducing ndpost
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) {
        message('The number of cores requested, ', mc.cores,
                       ',\n exceeds the number of cores detected via detectCores() ',
                       'reducing to ', mc.cores.detected)
        mc.cores <- mc.cores.detected
    }

    mc.ndpost <- ceiling(ndpost/mc.cores)

    ##if(length(keeptestfits)==0) keeptestfits <- (length(x.test)>0)

    post.list <- list()
    ##print(paste0('Parallel BART with mc.cores=', mc.cores))
    for(h in 1:H) {
        for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
              crisk3.bart(x.train=x.train, y.train=y.train,
                         x.train2=x.train2, y.train2=y.train2,
                         x.train3=x.train3, y.train3=y.train3,
                         times=times, delta=delta, K=K, events=events,
                         x.test=x.test, x.test2=x.test2, x.test3=x.test3,
                         cond=cond, cond2=cond2,
                         sparse=sparse, theta=theta, omega=omega,
                         a=a, b=b, augment=augment,
                         rho=rho, rho2=rho2, rho3=rho3,
                         xinfo=xinfo, xinfo2=xinfo2, xinfo3=xinfo3,
                         usequants=usequants,
                         rm.const=rm.const, type=type,
                         k=k, power=power, base=base,
                         offset=offset, offset2=offset2, offset3=offset3,
                         tau.num=tau.num, ntree=ntree, numcut=numcut,
                         ndpost=mc.ndpost, nskip=nskip, id=id,
                         keepevery = keepevery, printevery=printevery)},
                         ## keeptestfits = keeptestfits,
              silent=(i!=1))
              ## to avoid duplication of output
              ## capture stdout from first posterior only
        }

        post.list[[h]] <- parallel::mccollect()
    }

    if((H==1 & mc.cores==1) | attr(post.list[[1]][[1]], 'class')!='crisk3bart')
        return(post.list[[1]][[1]])
    else {
        p <- ncol(post.list[[1]][[1]]$tx.train)

        ndpost. <- 0
        for(h in 1:H) {
            mc.cores. = length(post.list[[h]])
            if(mc.cores!=mc.cores.)
                warning(paste0(
                    'The number of items returned by mccollect is ', mc.cores.))
            ndpost. <- ndpost.+mc.cores.*mc.ndpost

            for(i in mc.cores.:1) {
            if(h==1 & i==mc.cores.) {
                post <- post.list[[1]][[mc.cores.]]
                ##post$ndpost <- H*mc.cores*mc.ndpost
                post$ndpost <- ndpost.

                old.text <- paste0(as.character(mc.ndpost), ' ',
                                   as.character(ntree), ' ', as.character(p))
                old.stop <- nchar(old.text)

                post$treedraws$trees <- sub(old.text,
                                            paste0(as.character(post$ndpost),
                                                   ' ', as.character(ntree),
                                                   ' ', as.character(p)),
                                            post$treedraws$trees)

                old.text <- paste0(as.character(mc.ndpost), ' ',
                                   as.character(ntree), ' ', as.character(p))
                old.stop2 <- nchar(old.text)

                post$treedraws2$trees <- sub(old.text,
                                            paste0(as.character(post$ndpost),
                                                   ' ', as.character(ntree),
                                                   ' ', as.character(p)),
                                            post$treedraws2$trees)

                old.text <- paste0(as.character(mc.ndpost), ' ',
                                   as.character(ntree), ' ', as.character(p))
                old.stop3 <- nchar(old.text)

                post$treedraws3$trees <- sub(old.text,
                                            paste0(as.character(post$ndpost),
                                                   ' ', as.character(ntree),
                                                   ' ', as.character(p)),
                                            post$treedraws3$trees)

            }
            else {
                if(keeptestfits) {
                    post$yhat.test <- rbind(post$yhat.test,
                                            post.list[[h]][[i]]$yhat.test)
                    post$yhat.test2 <- rbind(post$yhat.test2,
                                             post.list[[h]][[i]]$yhat.test2)
                    post$yhat.test3 <- rbind(post$yhat.test3,
                                             post.list[[h]][[i]]$yhat.test3)
                    post$prob.test <- rbind(post$prob.test,
                                            post.list[[h]][[i]]$prob.test)
                    post$prob.test2 <- rbind(post$prob.test2,
                                             post.list[[h]][[i]]$prob.test2)
                    post$prob.test3 <- rbind(post$prob.test3,
                                             post.list[[h]][[i]]$prob.test3)
                    post$cif.test <- rbind(post$cif.test,
                                           post.list[[h]][[i]]$cif.test)
                    post$cif.test2 <- rbind(post$cif.test2,
                                            post.list[[h]][[i]]$cif.test2)
                    post$cif.test3 <- rbind(post$cif.test3,
                                            post.list[[h]][[i]]$cif.test3)
                    post$surv.test <- rbind(post$surv.test,
                                            post.list[[h]][[i]]$surv.test)
                }

                post$varcount <- rbind(post$varcount,
                                       post.list[[h]][[i]]$varcount)
                post$varcount2 <- rbind(post$varcount2,
                                        post.list[[h]][[i]]$varcount2)
                post$varcount3 <- rbind(post$varcount3,
                                        post.list[[h]][[i]]$varcount3)
                post$varprob <- rbind(post$varprob,
                                      post.list[[h]][[i]]$varprob)
                post$varprob2 <- rbind(post$varprob2,
                                       post.list[[h]][[i]]$varprob2)
                post$varprob3 <- rbind(post$varprob3,
                                       post.list[[h]][[i]]$varprob3)

                post$treedraws$trees <- paste0(post$treedraws$trees,
                                               substr(post.list[[h]][[i]]$treedraws$trees, old.stop+2,
                                                      nchar(post.list[[h]][[i]]$treedraws$trees)))

                post$treedraws2$trees <- paste0(post$treedraws2$trees,
                                               substr(post.list[[h]][[i]]$treedraws2$trees, old.stop2+2,
                                                      nchar(post.list[[h]][[i]]$treedraws2$trees)))

                post$treedraws3$trees <- paste0(post$treedraws3$trees,
                                               substr(post.list[[h]][[i]]$treedraws3$trees, old.stop3+2,
                                                      nchar(post.list[[h]][[i]]$treedraws3$trees)))
            }

            post.list[[h]][[i]] <- NULL
            }

        if(keeptestfits) {
            post$prob.test.mean  <- apply(post$prob.test,  2, mean)
            post$prob.test2.mean <- apply(post$prob.test2, 2, mean)
            post$prob.test3.mean <- apply(post$prob.test3, 2, mean)
            post$cif.test.mean  <- apply(post$cif.test,  2, mean)
            post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
            post$cif.test3.mean <- apply(post$cif.test3, 2, mean)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)
        }

        post$varcount.mean  <- apply(post$varcount,  2, mean)
        post$varcount2.mean <- apply(post$varcount2, 2, mean)
        post$varcount3.mean <- apply(post$varcount3, 2, mean)
        post$varprob.mean  <- apply(post$varprob,  2, mean)
        post$varprob2.mean <- apply(post$varprob2, 2, mean)
        post$varprob3.mean <- apply(post$varprob3, 2, mean)

        attr(post, 'class') <- 'crisk3bart'

        return(post)
    }
}
