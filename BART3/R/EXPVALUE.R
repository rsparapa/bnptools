
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2024-2025 Robert McCulloch and Rodney Sparapani
## EXPVALUE.R

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


## version using an array object
## these are wasteful since there are many zeros
## surprisingly, testing shows this is more efficient than sparse arrays
EXPVALUE = function(trees,
                    x.test,
                    S,
                    x.train)
                    ##call=FALSE)## default to R vs. C++ code
{
    P = ncol(x.test)
    N <- nrow(x.train)
    ## mask = integer(P)
    ## for(i in 1:P) mask[i] = (i %in% S)*1
    H = nrow(x.test)
    X.test <- x.test
    if(any(is.na(X.test)))  
        for(h in 1:H)  ## settings
            for(i in S) { ## subset of interest
                value <- X.test[h, i]
                while(is.na(value)) {
                    value <- x.train[sample.int(N, 1), i]
                    X.test[h, i] <- value
                }
            }

    ## if(call) {
    ##     return(.Call('cEXPVALUE', trees, x.test, mask))
    ## } else {
        M = dim(trees)[1]
        T = dim(trees)[2]
        node.max = dim(trees)[3]
        ## M*T*node.max*5
        G = function(n) ## node
        {
            if(trees[i, j, n, 1]==2) return(trees[i, j, n, 4]) ## a leaf
            else { ## a branch
                v=trees[i, j, n, 2]
                c=trees[i, j, n, 3]
                n=2*n
                m=n+1
                if(v %in% S) {
                    if(X.test[h, v]<c) return(G(n))
                    else return(G(m))
                } else {
                    a=trees[i, j, n, 5]
                    b=trees[i, j, m, 5]
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
##    }
}

## ## version using a sparse array object
## ## supposed to be more efficient since there are many zeros
## ## however, testing shows otherwise
## EXPVALUE = function(trees, x.test, S)
## {
##     H = nrow(x.test)
##     M = dim(trees)[1]
##     T = dim(trees)[2]
##     G = function(n) ## node
##     {
##         k=trees[cbind(i, j, n, 1:4)]
##         if(k[1]==2) return(k[4]) ## a leaf
##         else { ## a branch
##             v=k[2]
##             c=k[3]
##             n=2*n
##             m=n+1
##             if(v %in% S) {
##                 if(x.test[h, v]<c) return(G(n))
##                 else return(G(m))
##             } else {
##                 a=trees[matrix(c(i, j, n, 5), nrow=1, ncol=4)]
##                 b=trees[matrix(c(i, j, m, 5), nrow=1, ncol=4)]
##                 return((a*G(n)+b*G(m))/(a+b))
##             }
##         }
##     }
##     A = matrix(nrow=M, ncol=H)
##     B = matrix(nrow=M, ncol=T)
##     for(h in 1:H) { ## settings
##         for(i in 1:M) ## samples
##             for(j in 1:T) ## trees
##                 B[i, j]=G(1)
##         A[ , h]=apply(B, 1, sum)
##     }
##     return(A)
## }

## version using a list of lists
## about 10% of the size of an array
## EXPVALUE = function(trees, x.test, S)
## {
##     H = nrow(x.test)
##     M = length(trees)
##     T = length(trees[[1]])
##     G = function(n) ## node
##     {
##         if(trees[[i]][[j]][[1]][n, 1]==2)
##             return(trees[[i]][[j]][[2]][n]) ## a leaf
##         else { ## a branch
##             v=trees[[i]][[j]][[1]][n, 2]
##             c=trees[[i]][[j]][[1]][n, 3]
##             n=2*n
##             m=n+1
##             if(v %in% S) {
##                 if(x.test[h, v]<c) return(G(n))
##                 else return(G(m))
##             } else {
##                 a=trees[[i]][[j]][[1]][n, 4]
##                 b=trees[[i]][[j]][[1]][m, 4]
##                 return((a*G(n)+b*G(m))/(a+b))
##             }
##         }
##     }
##     A = matrix(nrow=M, ncol=H)
##     B = matrix(nrow=M, ncol=T)
##     for(h in 1:H) { ## settings
##         for(i in 1:M) ## samples
##             for(j in 1:T) ## trees
##                 B[i, j]=G(1)
##         A[ , h]=apply(B, 1, sum)
##     }
##     return(A)
## }
