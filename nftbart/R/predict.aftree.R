## Copyright (C) 2022 Rodney A. Sparapani

## This file is part of nftbart.
## predict.aftree.R

## nftbart is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.

## nftbart is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author contact information
## Rodney A. Sparapani: rsparapa@mcw.edu

predict.aftree = function(
                       ## data
                       object,
                       ## predictions
                       events=NULL,
                       FPD=FALSE,
                       hazard=FALSE,
                       probs=c(0.025, 0.975),
                       take.logs=TRUE,
                       seed=NULL,
                       ## default settings
                       ndpost=nrow(object$mix.prop),
                       nclust=ncol(object$mix.prop),
                       ## etc.
                       ...)
{
    if(length(object$m.test)==0)
        stop('x.test has to be provided to AFTrees()')

    L=ncol(object$m.test)

    if(length(events)>0) {
        events.matrix=(class(events)[1]=='matrix')

        if(events.matrix) {
            if(L!=nrow(events))
                stop(paste('the events matrix must have', L, 'rows'))
            K=ncol(events)
        } else { K=length(events) }

        if(any(is.na(c(events))) || K<=0)
            stop('events must be a vector, or matrix, of non-missing times')

        if(take.logs) events=log(events)
    } else {
        K=0
        FPD=FALSE
    }
    
    draw.logt=(length(seed)>0)
    
    N=ncol(object$m.train)
    H=L %/% N

    if(FPD) {
        if(L!=(N*H))
            stop("Friedman's partial dependence function: rows of x.test must be a multiple of x.train")
            
        if(events.matrix)
            stop("Friedman's partial dependence function: can't be used with a matrix of events")
    }
            
    res=list()

    if(draw.logt) {
        set.seed(seed)
        res$logt.test=matrix(0, nrow=ndpost, ncol=L)
        for(i in 1:L) {
            for(k in 1:nclust) {
                res$logt.test[ , i]=res$logt.test[ , i]+
                    object$mix.prop[ , k]*
                    rnorm(ndpost, object$locations[ , k]+
                             object$m.test[ , i],
                          object$sigma)
            }
        }
        res$logt.test.mean=apply(res$logt.test, 2, mean)
    }

if(K>0) {
    res$surv.test=matrix(0, nrow=ndpost, ncol=K*L)
    if(hazard)
        res$haz.test =matrix(0, nrow=ndpost, ncol=K*L)
    events.=events
    for(i in 1:L) {
        if(events.matrix) events.=events[i, ]
        for(m in 1:K) {
            l=(i-1)*K+m
            z=events.[m]
            t=exp(z)
            for(k in 1:nclust) {
                res$surv.test[ , l]=res$surv.test[ , l]+
                    object$mix.prop[ , k]*
                    pnorm(z, object$locations[ , k]+
                             object$m.test[ , i],
                          object$sigma, FALSE)

                if(hazard)
                    res$haz.test[ , l]=res$haz.test[ , l]+
                    object$mix.prop[ , k]*
                    dnorm(z, object$locations[ , k]+
                             object$m.test[ , i],
                          object$sigma)/(t*object$sigma*
                    pnorm(z, object$locations[ , k]+
                             object$m.test[ , i],
                          object$sigma, FALSE))
            }
            ## if(hazard) for(k in 1:nclust)
            ##         res$haz.test[ , l]=res$haz.test[ , l]+
            ##         (object$mix.prop[ , k]*
            ##         dnorm(z, object$locations[ , k]+
            ##                  object$m.test[ , i],
            ##               object$sigma)/
            ##         (t*object$sigma*res$surv.test[ , l]))
        }
    }

    if(FPD) {
        res$surv.fpd = matrix(0, nrow=ndpost, ncol=H*K)
        if(hazard)
            res$haz.fpd  = matrix(0, nrow=ndpost, ncol=H*K)

        for(i in 1:H) {
            h=seq(1, K*N, K)+(i-1)*K*N
            for(j in 1:K) {
                k=(i-1)*K+j
                res$surv.fpd[ , k]=apply(res$surv.test[ , h+j-1], 1, mean)
                if(hazard)
                    res$haz.fpd[ , k] =apply(res$haz.test[ , h+j-1], 1, mean)
            }
        }

        res$surv.test = NULL
        if(hazard) res$haz.test = NULL
        res$surv.fpd.mean =apply(res$surv.fpd, 2, mean)
        res$surv.fpd.lower=apply(res$surv.fpd, 2, quantile,
                                 probs=min(probs))
        res$surv.fpd.upper=apply(res$surv.fpd, 2, quantile,
                                 probs=max(probs))
        if(hazard) {
            res$haz.fpd.mean =apply(res$haz.fpd, 2, mean)
            res$haz.fpd.lower=apply(res$haz.fpd, 2, quantile,
                                    probs=min(probs))
            res$haz.fpd.upper=apply(res$haz.fpd, 2, quantile,
                                    probs=max(probs))
        }
    } else {
        res$surv.test.mean=apply(res$surv.test, 2, mean)
        res$surv.test.lower=apply(res$surv.test, 2, quantile, probs=min(probs))
        res$surv.test.upper=apply(res$surv.test, 2, quantile, probs=max(probs))
        if(hazard) {
            res$haz.test.mean=apply(res$haz.test, 2, mean)
            res$haz.test.lower=apply(res$haz.test, 2, quantile, probs=min(probs))
            res$haz.test.upper=apply(res$haz.test, 2, quantile, probs=max(probs))
        }
    }
}
    
    return(res)
}
