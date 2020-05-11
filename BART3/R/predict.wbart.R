
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Robert McCulloch and Rodney Sparapani

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

predict.wbart <- function(object, newdata, mc.cores=1,
                          openmp=(mc.cores.openmp()>0),
                          mult.impute=4, seed=99, ...)
{
    ##if(class(newdata) != "matrix") stop("newdata must be a matrix")

    p <- length(object$treedraws$cutpoints)

    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    if(.Platform$OS.type != "unix" || openmp || mc.cores==1) call <- pwbart
    else call <- mc.pwbart

    if(length(object$mu)==0) object$mu=object$offset

    ##return(call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$mu, ...))

    miss=apply(is.na(newdata), 2, sum)
    names(miss)=names(object$treedraws$cutpoints)
    miss.=(sum(miss)>0)

    if(!miss.) mult.impute=1
    else {
        set.seed(seed)
        newdata.=newdata
        n=nrow(newdata)
    }

    pred<-list()

    ## hot deck missing imputation
    for(k in 1:mult.impute) {
        if(miss.) {
            if(k==1)
                warning(paste0('missing elements of x imputed with hot decking ',
                               mult.impute, ' times'))
            else newdata=newdata.

            for(i in 1:n)
                for(j in 1:p) {
                    if(is.na(newdata[i, j]) && miss[j]==n) {
                        h=length(object$treedraws$cutpoints[[j]])
                        if(h==1) newdata[i, j]=object$treedraws$cutpoints[[j]][1]
                        else newdata[i, j]=object$treedraws$cutpoints[[j]][sample.int(h, 1)]
                    }
                    else while(is.na(newdata[i, j])) {
                        h=sample.int(n, 1)
                        newdata[i, j]=newdata[h, j]
                    }
                }
        }

        pred[[k]]=call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$mu, ...)

        if(mult.impute>1) {
            if(k==1) {
                pred[[1]]=pred[[1]]/mult.impute
            } else {
                pred[[1]]=pred[[1]]+pred[[k]]/mult.impute
            }
        }
    }

    ## pred[[1]]$mult.impute=mult.impute
    ## pred[[1]]$miss=miss

    ## if(mult.impute>1) {
    ##     pred[[1]]$prob.test.mean = apply(pred[[1]]$prob.test, 2, mean)
    ##     pred[[1]]$surv.test.mean = apply(pred[[1]]$surv.test, 2, mean)
    ##     pred[[1]]$tx.test=tx.test
    ## }

    return(pred[[1]])
}
