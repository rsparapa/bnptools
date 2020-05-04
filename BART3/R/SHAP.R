
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

SHAP=function(treedraws, ## treedraws item returned from BART
              x.train,   ## x.train to estimate coverage
              x.test,    ## settings of x.test: only x.test[ , S] are used,
                         ## but they must all be given
              S)         ## indices of subset
## Shapley additive explanation (SHAP) partial dependence function
{
    for(v in 1:ncol(x.train))
        if(any(is.na(x.train[ , v])))
            stop(paste0('x.train column with missing values:', v))

    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values:', v))

    tc <- textConnection(treedraws$tree)
    trees <- read.table(file=tc, fill=TRUE,
                    row.names=NULL, header=FALSE,
                    col.names=c('node', 'var', 'cut', 'leaf'))
    close(tc)
    trees.=trees[1, ] ## the first row are constants
    M=trees.$node ## number of samples
    T=trees.$var  ## number of trees
    trees=trees[-1, ]
    H=nrow(trees) ## length of the data frame
    node.max=max(trees$node)
    tier.max=floor(log2(node.max))
    Trees=array(0, dim=c(M, T, node.max, 5))
    ## last index of Trees
    ## 1 for nodes: 1 is a branch and 2 is a leaf
    ## 2 for variable: R index (add 1 to C/C++ index)
    ## 3 for cut-point: R index
    ## 4 for leaf value
    ## 5 for the number of training subjects passing through this node
    ## (data science/machine learning often calls this coverage)
    Cover=matrix(TRUE, nrow=N, ncol=node.max)
    i=1 ## index of samples
    j=0 ## index of trees
    for(l in 1:H) { ## index of the data.frame
        if(is.na(trees$var[l])) { ## the header for this tree
            j=j+1
            C=trees$node[l] ## number of nodes
            B=0             ## number of branches
            L=0             ## number of leaves
            for(h in l+(1:C)) { ## rows for this node
                k=trees$node[h]
                Trees[i, j, k, 1]=2 ## default to a leaf
                L=L+1
                if((k%%2)==0) {
                    Trees[i, j, k/2, 1]=1 ## but it was really a branch
                    B=B+1
                    L=L-1
                }
            }
            if(C!=B+L)
                stop(paste0('Sample:', i, ', Tree:', j, ', Row:', l,
                            ', C:', C, ', B:', B, ', L:', L))
            for(h in l+(1:C)) {
                k=trees$node[h]
                if(Trees[i, j, k, 1]==1) { ## a branch
                    v=trees$var[h]+1
                    Trees[i, j, k, 2]=v ## variable
                    c=treedraws$cutpoints[[v]][trees$cut[h]+1]
                    Trees[i, j, k, 3]=c ## cut-point
                } else if(Trees[i, j, k, 1]==2) { ## a leaf
                    Trees[i, j, k, 4]=trees$leaf[h]
                }
            }
            for(k in sort(trees$node[l+(1:C)])) {
                if(Trees[i, j, k, 1]==1) { ## a branch
                    v=Trees[i, j, k, 2]
                    c=Trees[i, j, k, 3]
                    if(k==1) {
                        Cover[ , 2]= x.train[ , v]<c
                        Cover[ , 3]= x.train[ , v]>=c
                    } else {
                        Cover[ , 2*k]  = Cover[ , k] & x.train[ , v]<c
                        Cover[ , 2*k+1]= Cover[ , k] & x.train[ , v]>=c
                    }
                }
                Trees[i, j, k, 5]=sum(Cover[ , k])
            }
            if((j%%T)==0) {
                i=i+1
                j=0
            }
        }
    }

    EXPVALUE = function()
    {
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
        H = nrow(x.test)
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
