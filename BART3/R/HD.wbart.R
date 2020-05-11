
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

## Hot deck (HD) partial dependence function
 HD.wbart=function(object,   ## object returned from BART
                   x.train, ## x.train to estimate coverage
                   x.test,  ## settings of x.test: only x.test[ , S]
                            ## are used but they must all be given
                   S,       ## indices of subset
                   mc.cores=1,
                   seed=99)
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

    N=nrow(x.train)
    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, S]==x.test[j, S]))
                stop(paste0('Row ', i, ' and ', j,
                            ' of x.test are equal with respect to S'))

    miss=apply(is.na(x.train), 2, sum)
    names(miss)=names(object$treedraws$cutpoints)
    miss.=(sum(miss)>0)

    pred=hotdeck(x.train, x.test, S, object$treedraws, object$offset,
                 mc.cores=mc.cores)

    ## set.seed(seed)
    ## for(i in 1:Q) {
    ##     X.test = x.train[sample.int(N), ] ## hot deck/permute x.train
    ##     for(j in 1:P) {
    ##         if(!(j %in% S) & miss[j]>0) {
    ##             for(k in 1:N)
    ##                 if(is.na(X.test[k, j]))
    ##                     while(is.na(X.test[k, j]))
    ##                         X.test[k, j]=X.test[sample.int(N, 1), j]
    ##         } else if(j %in% S) {
    ##             X.test[ , j]=x.test[i, j]
    ##         }
    ##     }
    ##     pred.=apply(predict(object, X.test, mc.cores=mc.cores,
    ##                         mult.impute=1), 1, mean)
    ##     if(i==1)
    ##         pred=cbind(pred.)
    ##     else
    ##         pred=cbind(pred, pred.)
    ## }

    return(pred)
}

