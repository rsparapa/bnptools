
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2023 Robert McCulloch and Rodney Sparapani
## FPD.wbart.R

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
FPD.wbart=function(object,  ## object returned from BART
                   x.test,  ## settings of x.test
                   S,       ## indices of subset
                   x.train=object$x.train,
                   probs=c(0.025, 0.975),
                   mc.cores=getOption('mc.cores', 1L),
                   mult.impute=4L,
                   seed=99L)
{
    P = ncol(x.train)

    if(!all(S %in% 1:P))
        stop('some elements of S are not columns of x.train')

    L=length(S)
    x.test=cbind(x.test)
    if(P==ncol(x.test)) x.test=cbind(x.test[ , S])
    else if(L!=ncol(x.test)) 
        stop('length of S and number of columns in x.test are not equal')

    if(any(is.na(c(x.test)))) stop('x.test cannot contain missing values')

    if(P!=length(object$treedraws$cutpoints))
        stop(paste0('the number of columns in x.train and\n',
                    'the length of cutpoints are not the same'))

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, ]==x.test[j, ]))
                stop(paste0('Row ', i, ' and ', j, ' of x.test are equal'))

    set.seed(seed)
    X.test = x.train
    for(i in 1:Q) {
        for(j in 1:L) X.test[ , S[j]]=x.test[i, j]
        ## for(j in S)
        ##     X.test[ , j]=x.test[i, j]
        pred.=apply(predict(object, X.test, mc.cores=mc.cores,
                            mult.impute=mult.impute, seed=NA), 1, mean)
        if(i==1)
            pred=cbind(pred.)
        else
            pred=cbind(pred, pred.)
    }

    return(pred)
}

