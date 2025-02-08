
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2024-2025 Robert McCulloch and Rodney Sparapani
## NN.wbart.R

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

## Nearest Neighbors marginal (NN) function
NN.wbart=function(object,  ## object returned from BART
                   x.test,  ## settings of x.test
                   S,       ## indices of subset
                   nearest,
                   x.train=object$x.train,
                   probs=c(0.025, 0.975),
                   mc.cores=getOption('mc.cores', 1L),
                   mult.impute=5L,
                   seed=99L)
{
    P = ncol(x.train)
    L=length(S)

    if(class(S)[1] == 'character') {
        class. <- class(x.train)[1]
        if(class. == 'matrix') names. <- dimnames(x.train)[[2]]
        else if(class. == 'data.frame') names. <- names(x.train)
        S. <- 0
        for(i in 1:L) {
            if(S[i] %in% names.) S.[i] <- which(S[i] == names.)
            else stop(paste0(S[i], 
                          ' has NOT been found to be a column name of x.train'))
        }
        S <- S.
    }

    if(!all(S %in% 1:P))
        stop('some elements of S are not columns of x.train')

    if(class(nearest)[1] == 'character') {
        class. <- class(x.train)[1]
        if(class. == 'matrix') names. <- dimnames(x.train)[[2]]
        else if(class. == 'data.frame') names. <- names(x.train)
        nearest. <- 0
        for(i in 1:length(nearest)) {
            if(nearest[i] %in% names.) 
                nearest.[i] <- which(nearest[i] == names.)
            else stop(paste0(nearest[i], 
                          ' has NOT been found to be a column name of x.train'))
        }
        nearest <- nearest.
    }

    if(!all(nearest %in% 1:P))
        stop('some elements of nearest are not columns of x.train')

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
    N <- nrow(x.train)
    dummy <- logical(P)
    for(j in 1:P) 
        dummy[j] <- (length(object$treedraws$cutpoints[[j]])==1)
    for(i in 1:Q) {
        X.test = x.train
        NN <- 1:N
        for(l in 1:L) {
            j <- S[l]
            if(j %in% nearest) {
                if(dummy[j]) { 
                    NN <- NN[x.train[NN, j] == x.test[i, l]]
                } else {
                    grid <- unique(sort(x.test[ , l]))
                    low <- -Inf
                    h <- which(grid == x.test[i, l])
                    if(h>1) low <- grid[h-1]
                    high <- Inf
                    if(h<length(grid)) high <- grid[h+1]
                    NN <- NN[low<x.train[NN, j] & x.train[NN, j]<high]
                }
            }
            X.test[ , j]=x.test[i, l] 
        }
        X.test <- X.test[NN, ]
        pred.=cbind(apply(predict(object, X.test, mc.cores=mc.cores,
                            mult.impute=mult.impute, seed=NA), 1, mean))
        
        if(i==1)
            pred=list(yhat.test      =pred.,
                      yhat.test.mean =mean(pred.),
                      yhat.test.lower=quantile(pred., probs=min(probs)),
                      yhat.test.upper=quantile(pred., probs=max(probs)))
        else {
            pred$yhat.test      =cbind(pred$yhat.test, pred.)
            pred$yhat.test.mean =c(pred$yhat.test.mean, mean(pred.))
            pred$yhat.test.lower=c(pred$yhat.test.lower, quantile(pred., min(probs)))
            pred$yhat.test.upper=c(pred$yhat.test.upper, quantile(pred., max(probs)))
        }
        
        ## pred.=apply(predict(object, X.test, mc.cores=mc.cores,
        ##                     mult.impute=mult.impute, seed=NA), 1, mean)
        ## if(i==1)
        ##     pred=cbind(pred.)
        ## else
        ##     pred=cbind(pred, pred.)
    }

    return(pred)
}

