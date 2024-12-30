
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2024 Robert McCulloch and Rodney Sparapani
## NN.survbart.R

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
NN.survbart=function(object,  ## object returned from BART
                   x.test,  ## settings of x.test
                   S,       ## indices of subset
                   subset.=NULL,
                   x.train=object$tx.test,
                   probs=c(0.025, 0.975),
                   mc.cores=getOption('mc.cores', 1L),
                   mult.impute=4L,
                   seed=99L)
{
    ## for(v in S)
    ##     if(any(is.na(x.test[ , v])))
    ##         stop(paste0('x.test column with missing values:', v))

    P = ncol(x.train)
    L=length(S)

    if(class(S)[1] == 'character') {
        class. <- class(x.train)[1]
        if(class. == 'matrix') names. <- dimnames(x.train)[[2]]
        else if(class. == 'data.frame') names. <- names(x.train)
        S. <- 0
        for(i in 1:L) {
            if(S[i] %in% names.) S.[i] <- which(S[i] == names.)
            else stop(paste0(S[i], ' has NOT been found to be a column name of x.train'))
        }
        S <- S.
    }
    
    if(!all(S %in% 1:P))
        stop('some elements of S are not columns of x.train')

    ## if(P!=ncol(x.test))
    ##     stop('the number of columns in x.train and x.test are not the same')
    ##L=length(S)
    x.test=cbind(x.test)
    if(P==ncol(x.test)) x.test=cbind(x.test[ , S])
    else if(L!=ncol(x.test)) 
        stop('length of S and number of columns in x.test are not equal')

    if(any(is.na(c(x.test)))) stop('x.test cannot contain missing values')

    if(P!=(length(object$treedraws$cutpoints)))
        stop(paste0('the number of columns in x.train and\n',
                    'the length of cutpoints are not the same'))

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, ]==x.test[j, ]))
                stop(paste0('Row ', i, ' and ', j, ' of x.test are equal'))

    NK=nrow(x.train)
    K=object$K
    N=NK/K
    M=object$ndpost
    set.seed(seed)
    ##X.test = x.train
    for(i in 1:Q) {
        X.test = x.train
        for(j in 1:L) {
            if(j %in% subset.) {
                ## assuming an increasing grid or a constant
                if(i==1) low=-Inf
                else low=x.test[i-1, j]
                if(low>x.test[i, j]) low=-Inf
                if(i==Q) high=Inf
                else high=x.test[i+1, j]
                if(high<x.test[i, j]) high=Inf
                if(low==x.test[i, j] | high==x.test[i, j])
                    X.test=X.test[X.test[ , S[j]]==x.test[i,j], ]
                else X.test=X.test[(low<X.test[ , S[j]] & X.test[ , S[j]]<high), ]
                ##print(c(low=low, high=high))
            }
            X.test[ , S[j]]=x.test[i, j] 
        }

        ##for(j in 1:L) X.test[ , S[j]]=x.test[i, j]
        ##for(j in S) X.test[ , j]=x.test[i, j]
        ## pre=surv.pre.bart(times=dots$times, delta=dots$delta,
        ##                   x.train=X.test, x.test=X.test,
        ##                   K=K, events=dots$events,
        ##                   ztimes=dots$ztimes, zdelta=dots$zdelta)
        ##pred=predict(object, pre$tx.test, mc.cores=mc.cores,
        pred=predict(object, X.test, mc.cores=mc.cores,
                            mult.impute=mult.impute, seed=NA)
        prob.test.=matrix(nrow=M, ncol=K)
        surv.test.=matrix(nrow=M, ncol=K)
        for(k in 1:K) {
            h=seq(k, NK, K)
            prob.test.[ , k]=apply(pred$prob.test[ , h], 1, mean)
            surv.test.[ , k]=apply(pred$surv.test[ , h], 1, mean)
        }
        if(i==1) {
            prob.test=cbind(prob.test.)
            surv.test=cbind(surv.test.)
        } else {
            prob.test=cbind(prob.test, prob.test.)
            surv.test=cbind(surv.test, surv.test.)
        }
    }
    prob.test.mean=apply(prob.test, 2, mean)
    prob.test.lower=apply(prob.test, 2, quantile, probs[1])
    prob.test.upper=apply(prob.test, 2, quantile, probs[2])
    surv.test.mean=apply(surv.test, 2, mean)
    surv.test.lower=apply(surv.test, 2, quantile, probs[1])
    surv.test.upper=apply(surv.test, 2, quantile, probs[2])
    pred=list(surv.test=surv.test, surv.test.mean=surv.test.mean,
              surv.test.lower=surv.test.lower, surv.test.upper=surv.test.upper,
              prob.test=prob.test, prob.test.mean=prob.test.mean,
              prob.test.lower=prob.test.lower, prob.test.upper=prob.test.upper)
    return(pred)
}

