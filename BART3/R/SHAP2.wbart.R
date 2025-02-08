
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2025 Robert McCulloch and Rodney Sparapani
## SHAP2.wbart.R

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
## for two-way interactions
SHAP2.wbart=function(object,  ## object returned from BART
              x.test,  ## settings of x.test
              S,       ## indices of two variables
              x.train=object$x.train,
              probs=c(0.025, 0.975),
              seed = 99L,
              mult.impute=5L
              )
{
    set.seed(seed)

    P = ncol(x.train)

    L=length(S)
    if(L!=2) stop('S must be of length 2')

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
        stop('some elements of S are not in x.train')

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
    
    ## if(P!=ncol(x.test))
    ##     stop('the number of columns in x.train and x.test are not the same')

    if(P!=length(object$treedraws$cutpoints))
        stop('the number of columns in x.train and length of cutpoints are not the same')

    Trees=read.trees(object$treedraws, x.train) 

    N <- nrow(x.train)

    i=S[1]
    j=S[2]
    
    M=P-1 ## use the same notation as SHAP
    L=M-1 ## but introduce L for convenience
    ## D=EXPVALUE(Trees, x.test, S, x.train)-
    ##     EXPVALUE(Trees, x.test, i, x.train)-
    ##     EXPVALUE(Trees, x.test, j, x.train)

    D <- 0

    for(k in 0:L) {
        for(h in 1:mult.impute) {
            x.test[ , -S] <- NA
            for(q in 1:Q) {
                missing <- is.na(x.test[q, ])
                while(any(missing)) {
                    x.test[q, missing] <- x.train[sample.int(N, 1), missing]
                    missing <- is.na(x.test[q, ])
                }
            }
                ## for(l in c(1:P)[-S]) { 
                ##     value <- x.test[q, l]
                ##     while(is.na(value)) {
                ##         value <- x.train[sample.int(N, 1), l]
                ##         x.test[q, l] <- value
                ##     }
                ## }
            if(h == 1) diff=0 
            if(k == 0) {
                diff <- diff+(EXPVALUE(Trees, x.test, S, x.train)-
                    EXPVALUE(Trees, x.test, i, x.train)-
                    EXPVALUE(Trees, x.test, j, x.train))/mult.impute
            } else {
                C=comb(L, k, (1:P)[-S])
                R=nrow(C)
                for(h in 1:R)
                    diff=diff+(EXPVALUE(Trees, x.test, c(C[h, ], S), x.train)-
                         EXPVALUE(Trees, x.test, c(C[h, ], i), x.train)-
                         EXPVALUE(Trees, x.test, c(C[h, ], j), x.train)+
                         EXPVALUE(Trees, x.test, C[h, ], x.train))/
                        (mult.impute*choose(L, k))
            }
        }
        D <- D+diff
    }

    return(0.5*D/M)
    
    ## unweighted for testing
    ## if(M>0) {
    ##     for(k in 1:M) {
    ##         C=comb(M, k, (1:P)[-S])
    ##         R=nrow(C)
    ##         N=N+R
    ##         for(h in 1:R)
    ##             D=D+EXPVALUE(Trees, x.test, c(C[h, ], S))-
    ##                 EXPVALUE(Trees, x.test, c(C[h, ], i))-
    ##                 EXPVALUE(Trees, x.test, c(C[h, ], j))+
    ##                 EXPVALUE(Trees, x.test, C[h, ])
    ##     }
    ##     D=D/N
    ## }

    ## return(D)
}
