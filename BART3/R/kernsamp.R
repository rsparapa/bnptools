
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2023 Robert McCulloch and Rodney Sparapani

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

kernsamp = function(
                   x.train,           ## matrix for training
                   x.test,	      ## matrix to predict at
                   S,                 ## x.test columns to condition on
                   ## others will be kernel samplinged
                   treedraws,	      ## $treedraws
                   mu=0,	      ## mean to add on
                   transposed=FALSE,
                   mult.impute=1L,
                   kern.var=FALSE,    ## return kernel sampling variance
                   alpha=0.05,        ## kernel sampling symmetric credible interval
                   probs=c(0.025, 0.975),
                                      ## equal-tail asymmetric credible interval
                   mc.cores=1L,       ## mc.kernsamp only
                   nice=19L           ## mc.kernsamp only
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

    if(kern.var && mult.impute==1)
        stop(paste0('To calculate kernel sampling variance, set mult.impute>1'))

    mask = 1:p
    for(i in mask)
        mask[i]=1*(i %in% S)

    res = as.list(1:mult.impute)
    for(i in 1:mult.impute) {
        res[[i]] = .Call("ckernsamp",
                    x.train,          ##training
                    x.test,           ##testing
                    as.integer(mask), ## 1 condition, 0 kernel sampling
                    treedraws         ##trees
                    )
        if(i==1) pred=res[[1]]
        else pred$yhat.test=pred$yhat.test+res[[i]]$yhat.test
    }
    pred$yhat.test=pred$yhat.test/mult.impute
    pred$yhat.test.mean=apply(pred$yhat.test, 2, mean)

    if(kern.var) {
        pred$yhat.test.var =apply(pred$yhat.test, 2, var)
        Yhat.test.mean=matrix(pred$yhat.test.mean, byrow=TRUE,
                              nrow=nrow(pred$yhat.test),
                              ncol=ncol(pred$yhat.test))
        for(i in 1:mult.impute) {
            if(i==1) pred$kern.test.var=(res[[1]]$yhat.test-Yhat.test.mean)^2
            else pred$kern.test.var=pred$kern.test.var+
                     ((res[[i]]$yhat.test-Yhat.test.mean)^2)
        }
        pred$kern.test.var=apply(pred$kern.test.var/(mult.impute^2), 2, mean)
        pred$yhat.test.var.=pred$yhat.test.var-pred$kern.test.var
        z=which(pred$yhat.test.var.<0)
        if(length(z)>0) {
            ## warning(paste0(length(z), ' adjusted kernel sampling variances<0 : ',
            ##                'increase mult.impute>', mult.impute))
            pred$yhat.test.var.[z]=pred$yhat.test.var[z]
        }
        pred$yhat.test=pred$yhat.test+mu
        pred$yhat.test.mean=pred$yhat.test.mean+mu
        ## z=qnorm(1-alpha/2)
        ## pred$yhat.test.lower.=pred$yhat.test.mean-z*sqrt(pred$yhat.test.var.)
        ## pred$yhat.test.upper.=pred$yhat.test.mean+z*sqrt(pred$yhat.test.var.)
        pred$yhat.test.=Yhat.test.mean+(pred$yhat.test-Yhat.test.mean)*
            sqrt(pred$yhat.test.var./pred$yhat.test.var)
        pred$yhat.test.lower=apply(pred$yhat.test., 2, quantile, probs=probs[1])
        pred$yhat.test.upper=apply(pred$yhat.test., 2, quantile, probs=probs[2])
        ##return(pred)
    } else {
        pred$yhat.test=pred$yhat.test+mu
        pred$yhat.test.mean=pred$yhat.test.mean+mu
        pred$yhat.test.lower=apply(pred$yhat.test, 2, quantile, probs=probs[1])
        pred$yhat.test.upper=apply(pred$yhat.test, 2, quantile, probs=probs[2])
        ##return(pred$yhat.test+mu)
    }
    return(pred)
}
