
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

hotdeck = function(
                   x.train,           ## matrix for training
                   x.test,	      ## matrix to predict at
                   S,                 ## subset of x.test columns to condition on
                                      ## others will be hot-decked
                   treedraws,	      ## $treedraws
                   mu=0,	      ## mean to add on
                   mc.cores=1L,       ## thread count
                   transposed=FALSE,
                   nice=19L           ## mc.hotdeck only
                   )
{
    if(!transposed) x.train <- t(bartModelMatrix(x.train))
    if(!transposed) x.test <- t(bartModelMatrix(x.test))

    p <- length(treedraws$cutpoints)

    if(p!=nrow(x.train))
        stop(paste0('The number of columns in x.train must be equal to ', p))

    if(p!=nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))

    mask = 1:p
    for(i in mask)
        mask[i]=1*(i %in% S)

    res = .Call("chotdeck",
                x.train,          ##training
                x.test,           ##testing
                as.integer(mask), ## 1 condition, 0 hot deck
                treedraws,        ##trees
                mc.cores          ##thread count
                )

    return(res$yhat.test+mu)
}
