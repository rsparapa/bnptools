
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
                     x.train=matrix(nrow=0, ncol=0),
                               ## x.train to estimate coverage
                     xinfo=NULL, numcut=100L)
{
    if(length(x.train)==0) stop('"x.train" must be specified')

    if(length(xinfo)==0) {
        bMM = bartModelMatrix(x.train, numcut=numcut)
        xinfo=bMM$xinfo
        numcut=bMM$numcut
    }

    N=nrow(x.train)
    P=ncol(x.train)
    if(numcut<P && length(numcut)==1) numcut=rep(numcut, P)
    if(length(numcut)!=P) stop(paste0('P=', P, ', numcut=', numcut))
    coverage=(N>0)
    if(coverage) {
        for(v in 1:ncol(x.train))
            if(any(is.na(x.train[ , v])))
                stop(paste0('x.train column with missing values:', v))
    }

    string=paste(1, ntree, paste0(P, '\n'))

    for(i in 1:ntree) {
        A=getTree(obj, i)
        h = nrow(A)
        string=paste0(string, paste(h, NA, NA, NA), '\n')
        if(h>0)
            for(j in 1:h) {
                if(A[j, 1]>0) { ## branch
                    v=A[j, 3]
                    c=A[j, 4]
                    ## RF is not using a grid: pick the nearest cutpoint
                    abs.diff=abs(xinfo[v, 1:numcut[v]]-c)
                    k=which(min(abs.diff)==abs.diff)
                    C=xinfo[v, k]
                    string=paste0(string, paste(j, v-1, k-1, 0), '\n')
                } else { ## leaf
                    string=paste0(string,
                                  paste(j, 0, 0,
                                        format(A[j, 6], digits=22)), '\n')
                }
            }
    }

    return(string)
}
