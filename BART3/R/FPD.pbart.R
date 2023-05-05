
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani
## FPD.pbart

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

## Friedman's partial dependence (FPD) function
FPD.pbart=function(object,  ## object returned from BART
##                 x.train, ## x.train to estimate coverage
                   x.test,  ## settings of x.test: only x.test[ , S]
                            ## are used but they must all be given
                   S,       ## indices of subset
                   x.train=object$x.train, ## x.train to estimate coverage
                   dots=NULL,## list of extra parameters if needed
                   mc.cores=getOption('mc.cores', 1L),
                   mult.impute=4L,
                   seed=99L)
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
        stop(paste0('the number of columns in x.train and\n',
                    'the length of cutpoints are not the same'))

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, S]==x.test[j, S]))
                stop(paste0('Row ', i, ' and ', j,
                            ' of x.test are equal with respect to S'))

    set.seed(seed)
    X.test = x.train
    for(i in 1:Q) {
        for(j in S)
            X.test[ , j]=x.test[i, j]
        pred=apply(predict(object, X.test, mc.cores=mc.cores,
                           mult.impute=mult.impute, seed=NA)$yhat.test,
                   1, mean)
        if(i==1) 
            yhat.test=cbind(pred)
        else
            yhat.test=cbind(yhat.test, pred)
    }

    prob.test=pnorm(yhat.test)
    yhat.test.mean=apply(yhat.test, 2, mean)
    prob.test.mean=apply(prob.test, 2, mean)
    pred=list(prob.test=prob.test, prob.test.mean=prob.test.mean,
              yhat.test=yhat.test, yhat.test.mean=yhat.test.mean)

    return(pred)
}

