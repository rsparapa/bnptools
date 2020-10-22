
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
                     x.train=matrix(nrow=0, ncol=0), xinfo=NULL)
{

    if(length(xinfo)==0) {
        if(length(x.train)==0) stop('"x.train" must be specified')
        bMM = bartModelMatrix(x.train, numcut=numcut)
        xinfo=bMM$xinfo
        numcut=bMM$numcut
        P=nrow(xinfo)
    } else {
        P=nrow(xinfo)
        numcut=integer(P)
        for(i in 1:P) numcut[i]=sum(!is.na(xinfo[i, ]))
    }

    if(length(numcut)!=P) stop(paste0('P=', P, ', numcut=', numcut))

    string=paste(1, ntree, paste0(P, '\n'))

    for(i in 1:ntree) {
        A=getTree(obj, i)
        h = nrow(A)
        string=paste0(string, h, '\n')
        if(h>0)
            for(j in 1:h) {
                if(A[j, 1]>0) { ## branch
                    v=A[j, 3]
                    if(numcut[v]==1) k=1
                    else {
                        ## RF is not using a grid: pick the nearest cutpoint
                        c=A[j, 4]
                        abs.diff=abs(xinfo[v, 1:numcut[v]]-c)
                        k=which(min(abs.diff)==abs.diff)
                    }
                    ##C=xinfo[v, k]
                    string=paste0(string, paste(j, v-1, k-1, 0), '\n')
                } else { ## leaf
                    string=paste0(string,
                                  paste(j, 0, 0,
                                        format(A[j, 6]/ntree, digits=22)), '\n')
                }
            }
    }

    return(string)
}
