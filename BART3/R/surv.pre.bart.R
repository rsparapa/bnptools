
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2025 Robert McCulloch and Rodney Sparapani

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


## you call this function before surv.bart()
## this function takes traditional time/delta
## survival variables and regressors (if any)
## and it constructs the corresponding
## tx.train, y.train and tx.test appropriate
## for use with pbart()

surv.pre.bart <- function(
                      times,
                      ## vector of survival times

                      delta,
                      ## vector of event indicators: 1 event, 0 censoring

                      x.train=NULL,
                      ## matrix of covariate regressors
                      ## can be NULL, i.e. KM analog

                      x.test=NULL,
                      ## matrix of covariate regressors at tx.test settings

                      K=NULL,
                      ## if specified, then use K quantiles for time grid

                      events=NULL,
                      ## if specified, then use events for time grid

                      ztimes=NULL,
                      zdelta=NULL,
                      ## column numbers of (ztimes, zdelta) time-dependent pairs

                      zsum = NULL,
                      ## list of time-dependent covariates to sum

                      rm.const=TRUE,
                      numcut=100,
                      grp=NULL, 
                      xinfo=matrix(0,0,0),
                      usequants=FALSE
                      ## parameters for bartModelMatrix
                      ) {
    ##binaryOffset <- qnorm(1-exp(-sum(delta)/sum(times)))

    if(length(numcut) == 1) numcut. <- numcut
    else numcut. <- 100

    N <- length(times)

    if(N!=length(delta))
        stop('The length of times and delta must be identical')

    if(length(x.train)>0 && N!=nrow(x.train))
        stop('The length of times and the number of rows in x.train, if any, must be identical')

    L <- length(ztimes)

    if(L!=length(zdelta))
        stop('The length of ztimes and zdelta, if any, must be identical')

    if(length(K)>0 || length(events)>0) {
        if(length(events)==0)
            events <- unique(quantile(times, probs=(1:K)/K))
        else if(!all(events==unique(events)))
            stop(paste0('events must be unique: ', events))
        attr(events, 'names') <- NULL

        events=sort(events)
        K <- length(events)

        for(i in 1:N) {
            if(times[i]>events[K]) {
                delta[i]=0
                times[i]=events[K]
            } else {
                k <- min(which(times[i]<=events))
                times[i] <- events[k]
            }
        }
    }
    else {
        events <- unique(sort(times))
        ## time grid of events including censoring times
        K <- length(events)
    }

    ##K <- length(events)

    if(events[1]<=0)
        stop('Time points exist less than or equal to time zero.')

    ## if(events[1]<0)
    ##     stop('Time points exist less than time zero.')
    ## else if(events[1]==0) {
    ##     warning('Time points exist equal to time zero.')
    ##     events=events[-1]
    ##     K=K-1
    ## }


    y.train <- integer(N) ## y.train is at least N long

    k <- 1

    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
        y.train[k] <- delta[i]*(times[i] == events[j])

        k <- k+1
    }

    m <- length(y.train)

    ##binaryOffset <- qnorm(mean(y.train))

    ## if(length(u.train)>0) {
    ##     makeU = TRUE
    ##     U.train <- integer(m)
    ## }
    ## else {
    ##     makeU = FALSE
    ##     U.train = NULL
    ## }

    if(length(x.train)==0) {
        p <- 0
        n <- 1

        X.train <- matrix(nrow=m, ncol=1, dimnames=list(NULL, 't'))
        Z <- 0
    } else {
        temp = bartModelMatrix(x.train, numcut=numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = (temp$X)
        nc = temp$numcut
        xinfo = temp$xinfo
        p <- ncol(x.train)
        
        if(length(x.test)>0) {
            x.test = (bartModelMatrix(x.test))
            if(class(rm.const)[1]=='logical' && rm.const)
                x.test = cbind(x.test[ , temp$rm.const])
            n <- nrow(x.test)
        }
        rm.const <- temp$rm.const
        if(length(grp)==0) grp <- temp$grp
        ##rm.const <- c(1, temp$rm.const+1)
        ##if(length(grp)==0) grp <- c(1, temp$grp)
        rm(temp)

        temp = bartModelMatrix(cbind(times), numcut=numcut, usequants=usequants)

        xinfo = rbind(temp$xinfo, xinfo)
        numcut = c(temp$numcut, nc)
        
        ## if(class(x.train)[1]=='data.frame') x.train=bartModelMatrix(x.train)

        ## p <- ncol(x.train)

        ## if(length(x.test)>0) {
        ##     if(class(x.test)[1]=='data.frame') x.test=bartModelMatrix(x.test)
        ##     n <- nrow(x.test)
        ## }

        Z <- length(zsum)

        X.train <- matrix(nrow=m, ncol=p+Z+1)

        znames <- names(zsum)
        if(Z>0 && length(znames) == 0) znames <- paste0('zsum', 1:Z)

        if(length(dimnames(x.train)[[2]])>0)
            dimnames(X.train)[[2]] <- c('times', dimnames(x.train)[[2]], znames)
        else dimnames(X.train)[[2]] <- c('times', paste0('x', 1:p), znames)
    }

    if(length(grp)==p) grp <- c(1, grp, rep(1, Z))
    if(length(grp)<p) grp <- rep(1, p+Z+1)
    if(length(rm.const) == p) rm.const <- c(1, rm.const+1)

    k <- 1
    
    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
        ##if(makeU) U.train[k] <- u.train[i]
        if(p==0) X.train[k, ] <- c(events[j])
        else X.train[k, ] <- c(events[j], x.train[i, ], rep(0, Z))

        k <- k+1
    }

    if(p==0 | length(x.test)>0) {
        X.test <- matrix(nrow=K*n, ncol=p+Z+1, dimnames=dimnames(X.train))

        for(i in 1:n) for(j in 1:K) {
            if(p==0) X.test[j, ] <- c(events[j])
            else X.test[(i-1)*K+j, ] <- c(events[j], x.test[i, ], rep(0, Z))
        }
    }
    else X.test <- matrix(nrow=0, ncol=0)*0

    if(L>0) {
        ztimes=ztimes+1
        zdelta=zdelta+1

        for(l in 1:L) {
            i=ztimes[l]
            j=zdelta[l]
            X.train[ , j]=(X.train[ , j]>0)*(X.train[ , 1]>=X.train[ , i])
            X.train[ , i]=X.train[ , 1]-(X.train[ , j]>0)*X.train[ , i]
            ##X.train[ , j]=X.train[ , j]*(X.train[ , 1]>=X.train[ , i])
            ##X.train[ , i]=X.train[ , 1]-X.train[ , j]*X.train[ , i]
            if(length(x.test)>0) {
                X.test[ , j]=(X.test[ , j]>0)*(X.test[ , 1]>=X.test[ , i])
                X.test[ , i]=X.test[ , 1]-(X.test[ , j]>0)*X.test[ , i]
                ##X.test[ , j]=X.test[ , j]*(X.test[ , 1]>=X.test[ , i])
                ##X.test[ , i]=X.test[ , 1]-X.test[ , j]*X.test[ , i]
            }
        }

        temp <- bartModelMatrix(cbind(X.train[,(p+2):(Z+p+1)]), numcut=numcut.)
        numcut <- c(numcut, temp$numcut)
        xinfo <- rbind(xinfo, temp$xinfo)
        rm.const <- c(rm.const, p+1+temp$rm.const)
    }

    if(Z>0) 
        for(l in 1:Z) {
            z <- 1+zsum[[l]]
            X.train[ , 1+p+l] <- apply(X.train[ , z], 1, sum)
            X.test[ , 1+p+l] <- apply(X.test[ , z], 1, sum)
        }

    return(list(y.train=y.train, tx.train=X.train, tx.test=X.test,
                times=events, K=K,
                rm.const=rm.const, grp=grp, xinfo=xinfo, numcut=numcut))
}
