
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

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

predict.crisk3bart <- function(object, newdata, newdata2, newdata3,
                               mc.cores=1,
                               openmp=(mc.cores.openmp()>0), ...) {

    p <- length(object$treedraws$cutpoints)
    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    p <- length(object$treedraws2$cutpoints)
    if(p!=ncol(newdata2))
        stop(paste0('The number of columns in newdata2 must be equal to ',p))

    p <- length(object$treedraws3$cutpoints)
    if(p!=ncol(newdata3))
        stop(paste0('The number of columns in newdata3 must be equal to ',p))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected)
        mc.cores <- mc.cores.detected

    if(length(object$binaryOffset)==0) object$binaryOffset=object$offset
    if(length(object$binaryOffset2)==0) object$binaryOffset2=object$offset2
    if(length(object$binaryOffset3)==0) object$binaryOffset3=object$offset3

    return(mc.crisk3.pwbart(newdata, newdata2, newdata3,
                           object$treedraws, object$treedraws2,
                           object$treedraws3,
                           object$binaryOffset, object$binaryOffset2,
                           object$binaryOffset3,
                           mc.cores=mc.cores, type=object$type, ...))
}

