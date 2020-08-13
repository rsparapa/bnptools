
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


mc.qbart <- function(x.train1=NULL, x.train2, times, delta,
                     x.test1=matrix(0,0,0), x.test2=matrix(0,0,0), K=100,
                     binoff=NA, grid=NA,
                     ntype=1,
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho1=NULL, rho2=NULL,
                     x1info=matrix(0,0,0), x2info=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     ##sigmaf=NA,
                     lambda=NA, ## tau.num=c(NA, 3, 6)[ntype],
                     ##tau.interval=0.9973,
                     ## offset=NULL,
                     w=rep(1, length(times)),
                     ntree1=50L, ntree2=200L, numcut1=100L, numcut2=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     mc.cores = 2L, nice = 19L, seed = 99L
                     )
{
    if(!(ntype %in% c(1,2))) stop('ntype must be either 1 or 2')

    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

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

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    mc.ndpost <- ceiling(ndpost/mc.cores)

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
            qbart(x.train1=x.train1, x.train2=x.train2, times=times, delta=delta,
                  x.test1=x.test1, x.test2=x.test2, K=K, binoff=binoff,grid=grid, ntype=ntype,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho1=rho1, rho2=rho2,
                  x1info=x1info, x2info=x2info, usequants=usequants,
                  rm.const=rm.const,
                  sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                  k=k, power=power, base=base,
                  ##sigmaf=sigmaf,
                  lambda=lambda, ## tau.num=tau.num,
                  ##tau.interval=tau.interval,
                  ## offset=offset,
                  w=w, ntree1=ntree1, ntree2=ntree2, numcut1=numcut1, numcut2=numcut2,
                  ndpost=mc.ndpost, nskip=nskip,
                  keepevery=keepevery, printevery=printevery,
                  transposed=TRUE)},
            silent=(i!=1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1 | attr(post, 'class')!="qbart") return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const <- rm.const

        post$ndpost <- mc.cores*mc.ndpost

        p1 <- nrow(x.train1[post$rm.const, ])

        old.text1 <- paste0(as.character(mc.ndpost), ' ', as.character(ntree1),
                           ' ', as.character(p1))
        old.stop1 <- nchar(old.text1)

        post$treedraws$trees1 <- sub(old.text1,
                                    paste0(as.character(post$ndpost), ' ',
                                           as.character(ntree1), ' ',
                                           as.character(p1)),
                                    post$treedraws$trees1)
        p2 <- nrow(x.train2[post$rm.const, ])

        old.text2 <- paste0(as.character(mc.ndpost), ' ', as.character(ntree2),
                           ' ', as.character(p2))
        old.stop2 <- nchar(old.text2)

        post$treedraws$trees2 <- sub(old.text2,
                                    paste0(as.character(post$ndpost), ' ',
                                           as.character(ntree2), ' ',
                                           as.character(p2)),
                                    post$treedraws$trees2)

        keeptest <- length(x.test1)>0

        for(i in 2:mc.cores) {
            post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$prob.train <- rbind(post$prob.train,
                                     post.list[[i]]$prob.train)

            post$y2hat.train <- rbind(post$y2hat.train,
                                     post.list[[i]]$y2hat.train)

            post$surv.train <- rbind(post$surv.train,
                                     post.list[[i]]$surv.train)
            if(keeptest) {
                post$prob.test <- rbind(post$prob.test,
                                     post.list[[i]]$prob.test)
                
                post$y2hat.test <- rbind(post$y2hat.test,
                                        post.list[[i]]$y2hat.test)

                post$surv.test <- rbind(post$surv.test,
                                        post.list[[i]]$surv.test)
            }

            post$varcount1 <- rbind(post$varcount1, post.list[[i]]$varcount1)
            post$varprob1 <- rbind(post$varprob1, post.list[[i]]$varprob1)
            post$varcount2 <- rbind(post$varcount2, post.list[[i]]$varcount2)
            post$varprob2 <- rbind(post$varprob2, post.list[[i]]$varprob2)

            post$treedraws$trees1 <- paste0(post$treedraws$trees1,
                                           substr(post.list[[i]]$treedraws$trees1, old.stop1+2,
                                                  nchar(post.list[[i]]$treedraws$trees1)))

            post$treedraws$trees2 <- paste0(post$treedraws$trees2,
                                           substr(post.list[[i]]$treedraws$trees2, old.stop2+2,
                                                  nchar(post.list[[i]]$treedraws$trees2)))

            post$proc.time['elapsed'] <- max(post$proc.time['elapsed'],
                                             post.list[[i]]$proc.time['elapsed'])
            for(j in 1:5)
                if(j!=3)
                    post$proc.time[j] <- post$proc.time[j]+post.list[[i]]$proc.time[j]
        }

        post$prob.train.mean <- apply(post$prob.train, 2, mean)
        post$y2hat.train.mean <- apply(post$y2hat.train, 2, mean)
        post$surv.train.mean <- apply(post$surv.train, 2, mean)

        if(keeptest) {
            post$prob.test.mean <- apply(post$prob.test, 2, mean)
            post$y2hat.test.mean <- apply(post$y2hat.test, 2, mean)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)
        }
        

        post$varcount1.mean <- apply(post$varcount1, 2, mean)
        post$varprob1.mean <- apply(post$varprob1, 2, mean)
        post$varcount2.mean <- apply(post$varcount2, 2, mean)
        post$varprob2.mean <- apply(post$varprob2, 2, mean)

        attr(post, 'class') <- 'qbart'

        return(post)
    }
}
