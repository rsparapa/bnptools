
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani
## write.trees.R

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

write.trees=function(treedraws, ## treedraws item returned from BART
                     thin = 1)  ## default to no thinning
{
    Trees. = read.trees(treedraws, call=TRUE)
    M = length(Trees.)              ## number of samples read in
    if(thin == M) seq.=M
    else seq. = seq(1, M, thin)
    M. = length(seq.) ## number of samples written out
    T = length(Trees.[[1]])         ## number of trees
    P = length(treedraws$cutpoints) ## number of variables
    line.=paste0(M., ' ', T, ' ', P, '\n')
    l = 1

    for(i in seq.)
        for(j in 1:T) {
            C = length(Trees.[[i]][[j]]$node)
            C. = length(which(Trees.[[i]][[j]]$node>0))
            l=l+1
            line.[l]=paste0(C., '\n')
            for(k in 1:C)
                if(Trees.[[i]][[j]]$node[k]==1) {
                    v = Trees.[[i]][[j]]$var[k]
                    c = Trees.[[i]][[j]]$cut[k]
                    l=l+1
                    line.[l]=paste0(k, ' ', v-1, ' ',
                                    which(c==treedraws$cutpoints[[v]])-1,
                                    ' -0.00000000000\n')
                } else if(Trees.[[i]][[j]]$node[k]==2) {
                    l=l+1
                    line.[l]=paste0(k, ' 0 0 ',
                                    Trees.[[i]][[j]]$leaf[k], '\n')
                }
        }
    ## if(file!="") cat(line., file=file, sep="")
    ## else return(line.)
    return(paste0(line., collapse=""))
}
