## Copyright (C) 2021 Rodney A. Sparapani

## This file is part of nftbart.
## bartModelMatrix.R

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

bartModelMatrix=function(X, numcut=0L, usequants=FALSE, type=7,
                         rm.const=FALSE, cont=FALSE, xicuts=NULL,
                         rm.vars=NULL) {
    if(length(rm.vars)>0) {
        ##rm.const=TRUE
        rm.vars[rm.vars>0]=-rm.vars[rm.vars>0]
    } else rm.vars=0

    X.class = class(X)[1]

    if(X.class=='factor') {
        X.class='data.frame'
        X=data.frame(X=X)
    }

    grp=0
    if(X.class=='data.frame') {
        p=dim(X)[2]
        xnm = names(X)
        j=0
        for(i in 1:p) {
            if(is.factor(X[[i]])) {
                Xtemp = class.ind(X[[i]])
                colnames(Xtemp) = paste(xnm[i],1:ncol(Xtemp),sep='')
                X[[i]]=Xtemp
                ##grp[i]=ncol(Xtemp)
                grp[i]=0
                for(k in j+(1:ncol(Xtemp)))
                    if(!(k %in% -rm.vars)) grp[i]=grp[i]+1
                j=j+ncol(Xtemp)
            } else {
                X[[i]]=cbind(X[[i]])
                colnames(X[[i]])=xnm[i]
                grp[i]=1-(j %in% -rm.vars)
                j=j+1
            }
        }
        names(grp)=xnm
        Xtemp=cbind(X[[1]])
        if(p>1) for(i in 2:p) Xtemp=cbind(Xtemp, X[[i]])
        X=Xtemp
    }
    else if(X.class=='numeric' | X.class=='integer') {
        X=cbind(as.numeric(X))
        grp=1
    }
    else if(X.class=='NULL') return(X)
    else if(X.class!='matrix')
        stop('Expecting either a factor, a vector, a matrix or a data.frame')

    N <- nrow(X)
    p <- ncol(X)

    ##xinfo. <- matrix(nrow=p, ncol=numcut)
    nc <- numcut
    ##rm.vars <- c()

    if(N>0 && p>0 && (rm.const || numcut[1]>0)) 
        for(j in 1:p) {
            X.class <- class(X[1, j])[1]

            if(X.class=='numeric' | X.class=='integer') {
                xs <- unique(sort(X[ , j]))
                k <- length(xs)
                nc[j] <- numcut

                if(cont) {
                    if(k==1) xs <-
                         seq(xs[1], xs[1]+1,
                             length.out=numcut+2)[-c(1, numcut+2)]
                    else xs <-
                         seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
                }
                else if(k==1) {
                    if(rm.const) rm.vars <- c(rm.vars, -j)
                    nc[j] <- 1
                }
                else if(k<numcut) {
                    xs <- 0.5*(xs[1:(k-1)]+xs[2:k])
                    nc[j] <- k-1
                }
                else if(usequants) {
                    xs <- quantile(X[ , j], type=type,
                                   probs=(0:(numcut+1))/
                                       (numcut+1))[-c(1, numcut+2)]
                    names(xs) <- NULL
                }
                else xs <-
                         seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
            }
            else
                stop(paste0('Variables of type ', X.class, ' are not supported'))

            ##nc[j] <- length(xs)
            ##xinfo.[j, 1:nc[j] ] <- xs
        }

    X <- data.matrix(X)

    if(numcut>0 && length(xicuts)==0)
        xicuts=xicuts(X, numcut=numcut)
    ## if(length(xinfo)>0) {
    ##     if(is.list(xinfo))
    ##         for(j in 1:p) xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
    ##     else if(is.matrix(xinfo)) xinfo. <- xinfo
    ##     else stop('Only a list or a matrix can be provided for xinfo')

    ##     for(j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
    ## }

    ##xinfo <- xinfo.

    if(rm.const || (length(rm.vars)>0 && rm.vars!=0)) {
        X <- X[ , rm.vars]
        nc <- nc[rm.vars]
        ##xinfo <- xinfo[rm.vars, ]
        for(i in length(rm.vars):1)
            xicuts[[-rm.vars[i] ]]=NULL
    }
    else if(length(rm.vars)==0) rm.vars <- 1:p

    ##dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)

    ## remove duplicate columns
    ## if(dedupe) {
    ##     na=is.na(apply(X, 1, sum))
    ##     XtX=crossprod(X[!na, ])
    ##     ev=eigen(XtX, TRUE, TRUE)$values
    ##     names(ev)=dimnames(X)[[2]]
    ## } else ev=NULL
    
    if(numcut==0) return(X)
    else {
        P=sum(grp)
        dummy=matrix(0, nrow=2, ncol=P)
        if(P==p) {
            h=1
            l=1
            for(i in 1:length(grp)) {
                for(j in 1:grp[i]) {
                    dummy[1, h]=l
                    dummy[2, h]=l+grp[i]-1
                    h=h+1
                }
                l=h
            }
            dimnames(dummy)[[2]]=dimnames(X)[[2]]
        }
        return(list(X=X, numcut=as.integer(nc), rm.const=rm.vars,
                    xicuts=xicuts, grp=grp, dummy=dummy))
    }
}
