
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
## for two-way interactions
SHAP2.wbart=function(object, ## object returned from BART
                    x.train, ## x.train to estimate coverage
                    x.test,  ## settings of x.test: only x.test[ , S]
                             ## are used but they must all be given
                    S)       ## indices of two variables
{
    if(length(S)!=2)
        stop('S must be of length 2')
    
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

    Trees=read.trees(object$treedraws, x.train)
    ## H = nrow(x.test)
    ## L = dim(Trees)[1]

    M=P-length(S)
    ## D=matrix(0, nrow=L, ncol=H)
    D=EXPVALUE(Trees, x.test, S)-
        EXPVALUE(Trees, x.test, i)-
        EXPVALUE(Trees, x.test, j) 
    N=1 ## number of terms
    i=S[1]
    j=S[2]
    
    if(M>0) {
        for(k in 1:M) {
            C=comb(M, k, (1:P)[-S])
            R=nrow(C)
            N=N+R
            for(h in 1:R)
                D=D+EXPVALUE(Trees, x.test, c(C[h, ], S))-
                    EXPVALUE(Trees, x.test, c(C[h, ], i))-
                    EXPVALUE(Trees, x.test, c(C[h, ], j))+
                    EXPVALUE(Trees, x.test, C[h, ])
        }
        D=D/N 
    }

    return(D)
}
