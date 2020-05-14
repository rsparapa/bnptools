
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


crisk3.bart <- function(
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
    ##keeptestfits = NULL,
    id = NULL,
    seed=99,    ## mc.crisk.bart only
    mc.cores=2, ## mc.crisk.bart only
    nice=19L    ## mc.crisk.bart only
)
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train3 <- bartModelMatrix(x.train3)
    x.test3 <- bartModelMatrix(x.test3)
    x.train2 <- bartModelMatrix(x.train2)
    x.test2 <- bartModelMatrix(x.test2)
    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    ##if(length(keeptestfits)==0) keeptestfits <- (length(x.test)>0)
    if(length(rho)==0) rho=ncol(x.train)
    if(length(rho2)==0) rho2=ncol(x.train2)
    if(length(rho3)==0) rho3=ncol(x.train3)

    if(length(y.train)==0) {
        pre <- crisk3.pre.bart(times, delta, x.train, x.test,
                              x.train2, x.test2, x.train3, x.test3,
                              K=K, events=events)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test
        y.train2 <- pre$y.train2
        x.train2 <- pre$tx.train2
        x.test2  <- pre$tx.test2
        y.train3 <- pre$y.train3
        x.train3 <- pre$tx.train3
        x.test3  <- pre$tx.test3

        times   <- pre$times
        K       <- pre$K

        if(length(cond)==0) cond <- pre$cond
        if(length(cond2)==0) cond2 <- pre$cond2
    }
    else {
        if(length(x.train)>0 & length(x.train2)>0 &
           nrow(x.train)!=nrow(x.train2))
            stop('number of rows in x.train and x.train2 must be equal')

        if(length(x.test)>0 & length(x.test2)>0 &
           nrow(x.test)!=nrow(x.test2))
            stop('number of rows in x.test and x.test2 must be equal')

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    keeptestfits <- (length(x.test)>0)

    if(length(xinfo)==0) {
        temp = bartModelMatrix(x.train3[cond2, ], numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train3 = t(temp$X)
        numcut3 = temp$numcut
        xinfo3 = temp$xinfo
        if(length(x.test3)>0) {
            x.test3 = t(bartModelMatrix(x.test3))
            if(class(rm.const)[1]=='logical' && rm.const)
                x.test3 = rbind(x.test3[temp$rm.const, ])
        }
        ##rm.const3 <- temp$rm.const
        rm(temp)

        temp = bartModelMatrix(x.train2[cond, ], numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train2 = t(temp$X)
        numcut2 = temp$numcut
        xinfo2 = temp$xinfo
        if(length(x.test2)>0) {
            x.test2 = t(bartModelMatrix(x.test2))
            if(class(rm.const)[1]=='logical' && rm.const)
                x.test2 = rbind(x.test2[temp$rm.const, ])
        }
        ##rm.const2 <- temp$rm.const
        rm(temp)

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

        xinfo2[1, ] <- xinfo[1, ] ## same time grid
        xinfo3[1, ] <- xinfo[1, ] ## same time grid
        transposed <- TRUE
    }
    else {
        x.train3=as.matrix(x.train3[cond2, ])
        x.train2=as.matrix(x.train2[cond, ])
        ##rm.const <- 1:ncol(x.train)
        ##rm.const2 <- 1:ncol(x.train2)
        ##rm.const3 <- 1:ncol(x.train3)
        transposed <- FALSE
    }

    post <- gbart(x.train=x.train, y.train=y.train,
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
                  transposed=transposed)
                  ##keeptestfits=keeptestfits

    if(type!=attr(post, 'class')) return(post)

    post2 <- gbart(x.train=x.train2, y.train=y.train2[cond],
                   x.test=x.test2, type=type,
                   sparse=sparse, theta=theta, omega=omega,
                   a=a, b=b, augment=augment, rho=rho2,
                   xinfo=xinfo2, usequants=usequants,
                   rm.const=rm.const,
                   k=k, power=power, base=base,
                   offset=offset2, tau.num=tau.num,
                   ntree=ntree, numcut=numcut2,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery,
                   printevery=printevery,
                   transposed=transposed)
                  ##keeptestfits=keeptestfits

    if(type!=attr(post2, 'class')) return(post2)

    post3 <- gbart(x.train=x.train3, y.train=y.train3[cond2],
                   x.test=x.test3, type=type,
                   sparse=sparse, theta=theta, omega=omega,
                   a=a, b=b, augment=augment, rho=rho3,
                   xinfo=xinfo3, usequants=usequants,
                   rm.const=rm.const,
                   k=k, power=power, base=base,
                   offset=offset3, tau.num=tau.num,
                   ntree=ntree, numcut=numcut3,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery,
                   printevery=printevery,
                   transposed=transposed)
                  ##keeptestfits=keeptestfits

    if(type!=attr(post3, 'class')) return(post3)

    post$offset2 <- post2$offset
    post$offset3 <- post3$offset
    post$id <- id
    post$times <- times
    post$K <- K

    if(!transposed) {
        post$tx.train <- x.train
        post$tx.train2 <- x.train2
        post$tx.train3 <- x.train3
    } else {
        post$tx.train <- t(x.train)
        post$tx.train2 <- t(x.train2)
        post$tx.train3 <- t(x.train3)
    }

    post$type <- type
    post$cond <- cond
    post$cond2 <- cond2
    post$treedraws2 <- post2$treedraws
    post$varcount2 <- post2$varcount
    post$varcount2.mean <- post2$varcount.mean
    post$varprob2 <- post2$varprob
    post$varprob2.mean <- post2$varprob.mean

    post$treedraws3 <- post3$treedraws
    post$varcount3 <- post3$varcount
    post$varcount3.mean <- post3$varcount.mean
    post$varprob3 <- post3$varprob
    post$varprob3.mean <- post3$varprob.mean

    post$rm.const <- rm.const
    ## post$rm.const2 <- rm.const2
    ## post$rm.const3 <- rm.const3
    post$yhat.train <- NULL
    post$yhat.train.mean <- NULL

    if(keeptestfits) {
        if(!transposed) {
            post$tx.test <- x.test
            post$tx.test2 <- x.test2
            post$tx.test3 <- x.test3
        } else {
            post$tx.test <- t(x.test)
            post$tx.test2 <- t(x.test2)
            post$tx.test3 <- t(x.test3)
        }

        H <- nrow(post$tx.test)/K ## the number of different settings

        post$yhat.test2 <- post2$yhat.test
        post$prob.test2 <- post2$prob.test
        post$prob.test2.mean <- post2$prob.test.mean

        post$yhat.test3 <- post3$yhat.test
        post$prob.test3 <- post3$prob.test
        post$prob.test3.mean <- post3$prob.test.mean

        post$surv.test <- (1-post$prob.test)*(1-post$prob.test2)*
            (1-post$prob.test3)
        post$cif.test <- post$prob.test
        post$cif.test2 <- (1-post$prob.test)*post$prob.test2
        post$cif.test3 <- (1-post$prob.test)*(1-post$prob.test2)*
            post$prob.test3

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$cif.test[ , l] <- post$cif.test[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test[ , l]
                post$cif.test2[ , l] <- post$cif.test2[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test2[ , l]
                post$cif.test3[ , l] <- post$cif.test3[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test3[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*
                    post$surv.test[ , l]
            }

        post$cif.test.mean <- apply(post$cif.test, 2, mean)
        post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
        post$cif.test3.mean <- apply(post$cif.test3, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'crisk3bart'

    return(post)
}
