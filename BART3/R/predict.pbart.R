
## BART: Bayesian Additive Regression Trees
## Copyright (C)-2020 2017 Robert McCulloch and Rodney Sparapani
## predict.pbart

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

predict.pbart <- function(object, newdata, mc.cores=1,
                          openmp=(mc.cores.openmp()>0),
                          mult.impute=4L, seed=99L, 
                          probs=c(0.025, 0.975), ...) {

    ##if(class(newdata) != "matrix") stop("newdata must be a matrix")

    p <- length(object$treedraws$cutpoints)

    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    if(.Platform$OS.type != "unix" || openmp || mc.cores==1) call <- pwbart
    else call <- mc.pwbart

    ##return(call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$binaryOffset, ...))

    if(length(object$binaryOffset)==0) object$binaryOffset=object$offset

    ## HD imputation 
    check = CDimpute(newdata)

    if(sum(check$miss.train)==0) mult.impute=1

    pred=list()
    
    for(i in 1:mult.impute) {
        yhat.test=call(newdata, object$treedraws, mc.cores=mc.cores,
                       mu=object$binaryOffset, ...)/mult.impute
        if(i==1) {
            sink('/dev/null')
            pred$yhat.test=yhat.test
        }
        else pred$yhat.test=pred$yhat.test+yhat.test
    }
    sink()
    
    pred$prob.test <- pnorm(pred$yhat.test)
    pred$prob.test.mean <- apply(pred$prob.test, 2, mean)
    pred$prob.test.lower <- apply(pred$prob.test, 2, quantile,
                                  probs=min(probs))
    pred$prob.test.upper <- apply(pred$prob.test, 2, quantile,
                                  probs=max(probs))
    pred$binaryOffset <- object$binaryOffset
    attr(pred, 'class') <- 'pbart'

    return(pred)
    
    ## pred <- list(yhat.test=call(newdata, object$treedraws, mc.cores=mc.cores,
    ##                             mu=object$binaryOffset, ...))

    ## pred$prob.test <- pnorm(pred$yhat.test)
    ## pred$prob.test.mean <- apply(pred$prob.test, 2, mean)
    ## pred$prob.test.lower <- apply(pred$prob.test, 2, quantile,
    ##                               probs=min(probs))
    ## pred$prob.test.upper <- apply(pred$prob.test, 2, quantile,
    ##                               probs=max(probs))
    ## pred$binaryOffset <- object$binaryOffset
    ## attr(pred, 'class') <- 'pbart'

    ## return(pred)
}

