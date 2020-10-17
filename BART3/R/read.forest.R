
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani
## read.forest.R

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

read.forest=function(obj, ## object returned from randomForest
                     ntree=200, maxnodes=4, node.max=2*maxnodes-1,
                     x.train=matrix(nrow=0, ncol=0))
                               ## x.train to estimate coverage
{
    N=nrow(x.train)
    P=ncol(x.train)
    coverage=(N>0)
    if(coverage) {
        for(v in 1:ncol(x.train))
            if(any(is.na(x.train[ , v])))
                stop(paste0('x.train column with missing values:', v))
    }

    string=paste(1, ntree, paste0(P, '\n'))
    tier.max=floor(log2(node.max))
    Trees=array(0, dim=c(ntree, node.max, 5))
    Trees. = list()
    ## last index of Trees
    ## 1 for nodes: 1 is a branch and 2 is a leaf
    ## 2 for variable: R index (add 1 to C/C++ index)
    ## 3 for cut-point
    ## 4 for leaf value
    ## 5 for the number of training subjects passing through this node
    ## (data science/machine learning often calls this coverage)
    if(coverage) Cover=matrix(TRUE, nrow=N, ncol=node.max)

    for(i in 1:ntree) {
        string=paste0(string, paste(node.max, NA, NA, NA), '\n')
        A=getTree(obj, i)
        h = nrow(A)
        Trees.[[i]]=list()
        Trees.[[i]]$node =integer(h)
        Trees.[[i]]$var  =integer(h)
        Trees.[[i]]$cut  =double(h)
        Trees.[[i]]$leaf =double(h)
        if(h>0)
            for(j in 1:h) {
                if(A[j, 1]>0) { ## branch
                    Trees[i, j, 1]=1
                    Trees.[[i]]$node[j]=1
                    Trees[i, j, 2]=A[j, 3]
                    Trees.[[i]]$var[j]=A[j, 3]
                    Trees[i, j, 3]=A[j, 4]
                    Trees.[[i]]$cut[j]=A[j, 4]
                    string=paste0(string, paste(j, A[j, 3]-1, A[j, 4], 0), '\n')
                } else { ## leaf
                    Trees[i, j, 1]=2
                    Trees.[[i]]$node[j]=2
                    Trees[i, j, 4]=A[j, 6]
                    Trees.[[i]]$leaf[j]=A[j, 6]
                    string=paste0(string, paste(j, 0, 0, A[j, 6]), '\n')
                }
            }
    }

    return(list(Trees=Trees, Trees.=Trees.))
}
