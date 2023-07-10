
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2023 Robert McCulloch and Rodney Sparapani
## SHAP.wbart.R

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
SHAP.wbart=function(object,  ## object returned from BART
                    x.test,  ## settings of x.test
              S,       ## indices of subset
              x.train=object$x.train,
              type='wbart',
              probs=c(0.025, 0.975),
              call=FALSE
              )
{
    ## for(v in S)
    ##     if(any(is.na(x.test[ , v])))
    ##         stop(paste0('x.test column with missing values:', v))

    P = ncol(x.train)

    if(!all(S %in% 1:P))
        stop('some elements of S are not in x.train')

    L=length(S)
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

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, ]==x.test[j, ]))
                warning(paste0('Row ', i, ' and ', j, ' of x.test are equal'))

    Trees=read.trees(object$treedraws, x.train, call)

    M=P-length(S)
    D=EXPVALUE(Trees, x.test, S, call) ## S vs. emptyset, i.e., subtract zero

    ## weighted difference
    if(M>0) {
        for(k in 1:M) {
            C=comb(M, k, (1:P)[-S])
            R=nrow(C)
            for(i in 1:R)
                D=D+(EXPVALUE(Trees, x.test, c(C[i, ], S), call)-
                    EXPVALUE(Trees, x.test, C[i, ], call))/choose(M, k)
        }
    }

    return(D/P)

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
