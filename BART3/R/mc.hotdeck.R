
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani

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

mc.hotdeck = function(
                   x.train,           ## matrix for training
                   x.test,	      ## matrix to predict at
                   S,                 ## x.test columns to condition on
                                      ## others will be hot-decked
                   treedraws,	      ## $treedraws
                   mu=0,	      ## mean to add on
                   transposed=FALSE,
                   mult.impute=1L,
                   hotd.var=FALSE,    ## hot-deck variance adjustment
                   alpha=0.05,        ## hot-deck symmetric credible interval
                   probs=c(0.025, 0.975),
                                      ## equal-tail asymmetric credible interval
                   mc.cores=2L,       ## mc.hotdeck only
                   nice=19L)          ## mc.hotdeck only
{
    ## if(.Platform$OS.type!='unix')
    ##     stop('parallel::mcparallel/mccollect do not exist on windows')

    if(!transposed) {
        x.train <- t(bartModelMatrix(x.train))
        x.test <- t(bartModelMatrix(x.test))
    }

    p <- length(treedraws$cutpoints)

    if(p!=nrow(x.train))
        stop(paste0('The number of columns in x.train must be equal to ', p))

    if(p!=nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))

    mc.cores.detected <- detectCores()

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected)
        mc.cores <- mc.cores.detected

    K <- ncol(x.test)
    k <- K%/%mc.cores-1
    j <- K
    for(i in 1:mc.cores) {
        if(i==mc.cores) h <- 1
        else h <- j-k

        parallel::mcparallel({psnice(value=nice);
            hotdeck(x.train=x.train,
                    x.test=matrix(x.test[ , h:j], nrow=p, ncol=j-h+1),
                    S=S, treedraws=treedraws, mu=0, transposed=TRUE,
                    mult.impute=mult.impute, hotd.var=hotd.var)},
            silent=(i!=1))
        j <- h-1
    }

    pred.list <- parallel::mccollect()

    if(hotd.var) {
        pred = list()
        pred$yhat.test <- pred.list[[1]]$yhat.test ## hot-deck posterior
        pred$yhat.test. <- pred.list[[1]]$yhat.test. ## adjusted posterior
        pred$hotd.test.var <- pred.list[[1]]$hotd.test.var

        if(mc.cores>1)
            for(i in 2:mc.cores) {
                pred$yhat.test <- cbind(pred$yhat.test,
                                        pred.list[[i]]$yhat.test)
                pred$yhat.test. <- cbind(pred$yhat.test.,
                                         pred.list[[i]]$yhat.test.)
                pred$hotd.test.var <- c(pred$hotd.test.var,
                                         pred.list[[i]]$hotd.test.var)
            }
        pred$yhat.test=pred$yhat.test+mu
        pred$yhat.test.=pred$yhat.test.+mu
        pred$yhat.test.mean=apply(pred$yhat.test, 2, mean)
        pred$yhat.test.var =apply(pred$yhat.test, 2, var)
        pred$yhat.test.var.=apply(pred$yhat.test., 2, var)
        z=qnorm(1-alpha/2)
        pred$yhat.test.lower.=pred$yhat.test.mean-z*sqrt(pred$yhat.test.var.)
        pred$yhat.test.upper.=pred$yhat.test.mean+z*sqrt(pred$yhat.test.var.)
        pred$yhat.test.lower=apply(pred$yhat.test., 2, quantile, probs=probs[1])
        pred$yhat.test.upper=apply(pred$yhat.test., 2, quantile, probs=probs[2])
        ##return(pred)
    } else {
        pred = list()
        pred$yhat.test <- pred.list[[1]]$yhat.test ## hot-deck posterior

        if(mc.cores>1)
            for(i in 2:mc.cores) {
                pred$yhat.test <- cbind(pred$yhat.test,
                                        pred.list[[i]]$yhat.test)
            }
        pred$yhat.test=pred$yhat.test+mu
        pred$yhat.test.mean=apply(pred$yhat.test, 2, mean)
        pred$yhat.test.lower=apply(pred$yhat.test, 2, quantile, probs=probs[1])
        pred$yhat.test.upper=apply(pred$yhat.test, 2, quantile, probs=probs[2])
        ## pred <- pred.list[[1]]

        ## if(mc.cores>1)
        ##     for(i in 2:mc.cores)
        ##         pred <- cbind(pred, pred.list[[i]])
        ## return(pred+mu)
    }
    return(pred)
}
