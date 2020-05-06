
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

## Shapley additive explanation (SHAP) partial dependence function
SHAP.survbart  =function(object,       ## object returned from BART
                         x.train,      ## x.train to estimate coverage
                         x.test,       ## settings of x.test: only x.test[ , S]
                                       ## are used but they must all be given
                         S,            ## indices of subset
                         type='pbart', ## type of probability model
                         probs=c(0.025, 0.975))
{
    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values:', v))

    P = ncol(x.train)

    if(!all(S %in% 1:P))
        stop('some elements of S are not in x.train')

    if(P!=ncol(x.test))
        stop('the number of columns in x.train and x.test are not the same')

    if(P!=length(object$treedraws$cutpoints))
        stop('the number of columns in x.train and length of cutpoints are not the same')

    trees=read.trees(object$treedraws, x.train)

    pred=list()

    pred$times <- object$times
    K <- object$K
    pred$K = K
    pred$yhat.test <- object$offset+EXPVALUE(trees, x.test, S)
    pred$offset <- object$offset

    if(type=='pbart')
        pred$prob.test <- pnorm(pred$yhat.test)
    else if(type=='lbart')
        pred$prob.test <- plogis(pred$yhat.test)

    pred$surv.test <- 1-pred$prob.test

    H <- nrow(x.test)/K ## the number of different settings

    for(h in 1:H)
        for(j in 2:K) {
            l <- K*(h-1)+j

            pred$surv.test[ , l] <- pred$surv.test[ , l-1]*
                pred$surv.test[ , l]
        }

    pred$surv.test.mean <- apply(pred$surv.test, 2, mean)
    pred$surv.test.lower <- apply(pred$surv.test, 2, quantile, probs=probs[1])
    pred$surv.test.upper <- apply(pred$surv.test, 2, quantile, probs=probs[2])

    attr(pred, 'class') <- 'survbart'

    return(pred)
}
