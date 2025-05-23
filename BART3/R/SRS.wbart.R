
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2024 Robert McCulloch and Rodney Sparapani
## SRS.wbart.R

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

## kernel sampling for Friedman's partial dependence (FPD) function
SRS.wbart=function(object, ## object returned from BART
                   x.test,  ## settings of x.test
                   S,       ## indices of subset
                   x.train=object$x.train, 
                   probs=c(0.025, 0.975),
                   mc.cores=getOption('mc.cores', 1L),
                   mult.impute=30L,
                   seed=99L,
                   kern.var=FALSE, ## kernel sampling variance adjustment
                   alpha=0.05, ## kernel sampling symmetric credible interval
                   nice=19L)
{
    if(mc.cores>1L) {
        if(.Platform$OS.type!='unix')
            stop('parallel::mcparallel/mccollect do not exist on windows')
        else {
            RNGkind("L'Ecuyer-CMRG")
            set.seed(seed)
            parallel::mc.reset.stream()
        }
    } else { set.seed(seed) }

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

    ##L=length(S)
    x.test=cbind(x.test)
    Q=nrow(x.test)
    if(L==ncol(x.test)) {
        X.test=x.test
        x.test=matrix(0, nrow=Q, ncol=P)
        for(j in 1:L) x.test[ , S[j]]=X.test[ , j]
    } else if(P!=ncol(x.test)) { 
        stop('the number of columns in x.train and x.test are not equal') }

    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values: S=', v))

    if(P!=length(object$treedraws$cutpoints))
        stop(paste0('the number of columns in x.train and\n',
                    'the length of cutpoints are not the same'))

    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, S]==x.test[j, S]))
                stop(paste0('Row ', i, ' and ', j,
                            ' of x.test are equal with respect to S'))

    if(mc.cores>1L)
        return(mc.kernsamp(x.train, x.test, S, object$treedraws,
                          object$offset, mult.impute=mult.impute,
                          kern.var=kern.var, alpha=alpha, probs=probs,
                          mc.cores=mc.cores, nice=nice))
    else
        return(kernsamp(x.train, x.test, S, object$treedraws,
                       object$offset, mult.impute=mult.impute,
                       kern.var=kern.var, alpha=alpha, probs=probs))
}

