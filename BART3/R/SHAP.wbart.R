
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
SHAP.wbart=function(object,  ## object returned from BART
                    x.train, ## x.train to estimate coverage
                    x.test,  ## settings of x.test: only x.test[ , S]
                             ## are used but they must all be given
                    S)       ## indices of subset
{
    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values:', v))

    Trees=read.trees(object$treedraws, x.train)

    EXPVALUE = function()
    {
        H = nrow(x.test)
        M = dim(Trees)[1]
        T = dim(Trees)[2]
        G = function(n) ## node
        {
            if(Trees[i, j, n, 1]==2) return(Trees[i, j, n, 4]) ## a leaf
            else { ## a branch
                v=Trees[i, j, n, 2]
                c=Trees[i, j, n, 3]
                n=2*n
                m=n+1
                if(v %in% S) {
                    if(x.test[h, v]<c) return(G(n))
                    else return(G(m))
                } else {
                    a=Trees[i, j, n, 5]
                    b=Trees[i, j, m, 5]
                    return((a*G(n)+b*G(m))/(a+b))
                }
            }
        }
        A = matrix(nrow=M, ncol=H)
        B = matrix(nrow=M, ncol=T)
        for(h in 1:H) { ## settings
            for(i in 1:M) ## samples
                for(j in 1:T) ## trees
                    B[i, j]=G(1)
            A[ , h]=apply(B, 1, sum)
        }
        return(A)
    }

    return(EXPVALUE())
}
