
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani
## read.trees.R

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

read.trees=function(treedraws, ## treedraws item returned from BART
                    x.train=matrix(nrow=0, ncol=0),
                               ## x.train to estimate coverage
                    call=FALSE,## default to R vs. C++ code
                    cutpoints=NULL,
                    trees=NULL)
{
    N=nrow(x.train)
    coverage=(N>0)
    if(coverage) {
        for(v in 1:ncol(x.train))
            if(any(is.na(x.train[ , v])))
                stop(paste0('x.train column with missing values:', v))
    }

    if(length(treedraws$cutpoints)==0) {
        if(length(cutpoints)>0)
            treedraws$cutpoints=cutpoints
        else stop('The cutpoints item was not found in treedraws')
    }
    if(length(treedraws$trees)==0) {
        if(length(trees)>0)
            treedraws$trees=trees
        else stop('The trees string was not found in treedraws')
    }
    
    ##print(paste0('tc <- textConnection(treedraws$', name., ')'))
    tc <- textConnection(treedraws$trees)
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
    if(call) Trees. = replicate(n=M, expr=list())
    ## last index of Trees
    ## 1 for nodes: 1 is a branch and 2 is a leaf
    ## 2 for variable: R index (add 1 to C/C++ index)
    ## 3 for cut-point
    ## 4 for leaf value
    ## 5 for the number of training subjects passing through this node
    ## (data science/machine learning often calls this coverage)
    if(coverage) Cover=matrix(TRUE, nrow=N, ncol=node.max)

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
                Trees[i, j, k, 1]=2L ## default to a leaf
                L=L+1
                if((k%%2)==0) {
                    Trees[i, j, k/2, 1]=1L ## but it was really a branch
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
            if(coverage) {
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
            }
            if(call) {
                node.max=max(which(Trees[i, j, , 1]==2L))
                Trees.[[i]][[j]]=list()
                Trees.[[i]][[j]]$node =integer(node.max)
                Trees.[[i]][[j]]$var  =integer(node.max)
                Trees.[[i]][[j]]$cut  =double(node.max)
                Trees.[[i]][[j]]$leaf =double(node.max)
                if(coverage) Trees.[[i]][[j]]$cover=integer(node.max)
                k=which(Trees[i, j, , 1]==1L)
                Trees.[[i]][[j]]$node[k] =1L ## branches
                Trees.[[i]][[j]]$var[k]  =Trees[i, j, k, 2] ## variables
                Trees.[[i]][[j]]$cut[k]  =Trees[i, j, k, 3] ## cut-points
                if(coverage) Trees.[[i]][[j]]$cover[k]=Trees[i, j, k, 5] ## coverage
                k=which(Trees[i, j, , 1]==2L) ## leaves
                Trees.[[i]][[j]]$node[k] =2L
                Trees.[[i]][[j]]$leaf[k] =Trees[i, j, k, 4] ## leaves
                if(coverage) Trees.[[i]][[j]]$cover[k]=Trees[i, j, k, 5] ## coverage
            }
            if((j%%T)==0) {
                i=i+1
                j=0
            }
        }
    }
    ##print(object.size(Trees))
    ##print(object.size(Trees.))
    if(call) return(Trees.)
    else return(Trees)
}

