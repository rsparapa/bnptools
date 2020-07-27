
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

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

hotdeck = function(
                   x.train,           ## matrix for training
                   x.test,	      ## matrix to predict at
                   S,                 ## x.test columns to condition on
                   ## others will be hot-decked
                   treedraws,	      ## $treedraws
                   mu=0,	      ## mean to add on
                   transposed=FALSE,
                   mult.impute=1L,
                   hdvar=FALSE,       ## return hot-deck variance
                   mc.cores=1L,       ## mc.hotdeck only
                   nice=19L           ## mc.hotdeck only
                   )
{
    if(!transposed) {
        x.train <- t(bartModelMatrix(x.train))
        x.test <- t(bartModelMatrix(x.test))
    }

    p <- length(treedraws$cutpoints)

    if(p!=nrow(x.train))
        stop(paste0('The number of columns in x.train must be equal to ', p))

    if(p!=nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))

    if(hdvar && mult.impute==1)
        stop(paste0('To calculate hot-deck variance, set mult.impute>1'))

    mask = 1:p
    for(i in mask)
        mask[i]=1*(i %in% S)

    res = as.list(1:mult.impute)
    for(i in 1:mult.impute) {
        res[[i]] = .Call("chotdeck",
                    x.train,          ##training
                    x.test,           ##testing
                    as.integer(mask), ## 1 condition, 0 hot deck
                    treedraws         ##trees
                    )
        if(i==1) pred=res[[1]]
        else pred$yhat.test=pred$yhat.test+res[[i]]$yhat.test
    }
    pred$yhat.test=pred$yhat.test/mult.impute
    pred$yhat.test.mean=apply(pred$yhat.test, 2, mean)
    yhat.test.mean=matrix(pred$yhat.test.mean, byrow=TRUE,
                          nrow=nrow(pred$yhat.test),
                          ncol=ncol(pred$yhat.test))
    if(hdvar) {
        for(i in 1:mult.impute) {
            if(i==1) pred$hdvar.test=(res[[1]]$yhat.test-yhat.test.mean)^2
            else pred$hdvar.test=pred$hdvar.test+
                     ((res[[i]]$yhat.test-yhat.test.mean)^2)
        }
        pred$hdvar.test=pred$hdvar.test/mult.impute
        pred$yhat.test=pred$yhat.test+mu
        pred$yhat.test.mean=pred$yhat.test.mean+mu
        return(pred)
    } else {
        return(pred$yhat.test+mu)
    }

    ## OpenMP will not work here since we cannot
    ## estimate 1, ..., ndpost predictions simultaneously
    ## for each sample, we need to hotdeck the values
    ## however, we can parallel-ize predictions 1, ..., np
    ## but we are not using OpenMP: just multiple threads
    ## each given a block of xtest which is even easier
    ## see mc.hotdeck
}
