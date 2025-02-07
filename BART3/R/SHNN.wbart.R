
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2025 Robert McCulloch and Rodney Sparapani
## SHNN.wbart.R

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

## Shapley additive explanation with Nearest Neighbors 
SHNN.wbart=function(object,  ## object returned from BART
                    x.test,  ## settings of x.test
              S=NULL,        ## indices of subset
              x.train=object$x.train,
              type='wbart',
              probs=c(0.025, 0.975),
              seed = 99L,
              mult.impute=5L
              )
{
    set.seed(seed)
    call <- FALSE

    P = ncol(x.train)
    L=length(S)

    ##if(L>2) stop('NN is only supports 1 or 2 variables of interest')

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

    if(P!=length(object$treedraws$cutpoints))
        stop('the number of columns in x.train and length of cutpoints are not the same')

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, ]==x.test[j, ]))
                warning(paste0('Row ', i, ' and ', j, ' of x.test are equal'))

    Trees=read.trees(object$treedraws, x.train, call)

    M=P-L

    if(M<=0) stop('The length of S must be smaller than P')

    ##dummy <- logical(L)
    ##for(l in 1:L) dummy[l] <- (length(object$treedraws$cutpoints[[ S[l] ]])==1)
    N <- nrow(x.train)
    dummy <- logical(P)
    for(j in 1:P) dummy[j] <- (length(object$treedraws$cutpoints[[j]])==1)

    P.=lfactorial(P)
    D <- 0
    ## D=EXPVALUE(Trees, x.test, S, call)*exp(lfactorial(M)-P.)
    ## S vs. emptyset, i.e., subtract zero

    ## weighted difference
    for(k in 0:M) { 
        for(j in 1:mult.impute) {
            x.test[ , -S] <- NA
            for(h in 1:Q) {
                missing <- is.na(x.test[h, ])
                while(any(missing)) {
                        NN <- 1:N
                        for(l in 1:L) {
                            if(dummy[S[l]]) { 
                               NN <- NN[x.train[NN, S[l]] == x.test[h, S[l]]]
                            } else {
                                low <- -Inf
                                high <- Inf
                                ## if(h>1 && x.test[h-1, S[l]]<x.test[h, S[l]])
                                ##     low <- x.test[h-1, S[l]]
                                ## if(h<Q && x.test[h, S[l]]<x.test[h+1, S[l]])
                                ##     high <- x.test[h+1, S[l]]
                                grid <- unique(sort(x.test[ , S[l]]))
                                i <- which(grid == x.test[h, S[l]])
                                if(i>1) low <- grid[i-1]
                                if(i<length(grid)) high <- grid[i+1]
                                NN <- NN[low<x.train[NN, S[l]] & x.train[NN, S[l]]<high]
                            }
                        }
                        x.test[h, missing] <- x.train[sample(NN, 1), missing]
                    missing <- is.na(x.test[h, ])
                    }
                }
            if(j == 1) diff=0 
            if(k == 0) {
                diff=diff+EXPVALUE(Trees, x.test, S, x.train)*
                    exp(lfactorial(M)-P.)/mult.impute 
            } else {
                C=comb(M, k, (1:P)[-S])
                R=nrow(C)
                for(i in 1:R)
                   diff <- diff+(EXPVALUE(Trees, x.test, c(C[i, ], S), x.train)-
                                  EXPVALUE(Trees, x.test, C[i, ], x.train))*
                        exp(lfactorial(k)+lfactorial(M-k)-P.)/mult.impute
            }
        }
        D <- D+diff
    }

    return(list(yhat.test=D,
           yhat.test.mean =apply(D, 2, mean),
           yhat.test.lower=apply(D, 2, quantile, probs=min(probs)),
           yhat.test.upper=apply(D, 2, quantile, probs=max(probs))))
    
    ## M=P-length(S)
    ## D=EXPVALUE(Trees, x.test, S, call) ## S vs. emptyset, i.e., subtract zero

    ## ## weighted difference
    ## if(M>0) {
    ##     for(k in 1:M) {
    ##         C=comb(M, k, (1:P)[-S])
    ##         R=nrow(C)
    ##         for(i in 1:R)
    ##             D=D+(EXPVALUE(Trees, x.test, c(C[i, ], S), call)-
    ##                 EXPVALUE(Trees, x.test, C[i, ], call))/choose(M, k)
    ##     }
    ## }

    ## D=D/P
    ## D=list(yhat.test=D,
    ##        yhat.test.mean =apply(D, 2, mean),
    ##        yhat.test.lower=apply(D, 2, quantile, probs=min(probs)),
    ##        yhat.test.upper=apply(D, 2, quantile, probs=max(probs)))
    ## return(D)

    ## unweighted difference
    ## N=1 ## number of terms

    ## if(M>0) {
    ##     for(k in 1:M) {
    ##         C=comb(M, k, (1:P)[-S])
    ##         R=nrow(C)
    ##         N=N+R
    ##         for(i in 1:R)
    ##             D=D+EXPVALUE(Trees, x.test, c(C[i, ], S))-
    ##                 EXPVALUE(Trees, x.test, C[i, ])
    ##     }
    ##     D=D/N
    ## }
}
