##hbartModelMatrix=function(xdf) {
## p=ncol(xdf)
## Xnum = NULL #this will be all the numeric columns
## Xfac = NULL #this will be all the factors expanded into dummies
## for(i in 1:p) { # for each xi add to Xnum if numeric else add dummies to Xfac
##    xnm = names(xdf)[i]
##    if(is.factor(xdf[[i]])) {
##       Xtemp = nnet::class.ind(xdf[[i]])
##       colnames(Xtemp) = paste(xnm,1:ncol(Xtemp),sep='')
##       Xfac = cbind(Xfac,Xtemp)
##    } else {
##       Xnum=cbind(Xnum,xdf[[i]])
##       colnames(Xnum)[ncol(Xnum)]=xnm
##    }
## }
## return(cbind(Xnum,Xfac))

hbartModelMatrix=function(xdf) {
    X=xdf
    numcut=0L
    usequants=FALSE
    type=7
    rm.const=FALSE
    cont=FALSE
    xinfo=NULL

    X.class = class(X)[1]

    if(X.class=='factor') {
        X.class='data.frame'
        X=data.frame(X=X)
    }

    grp=NULL

    if(X.class=='data.frame') {
        p=dim(X)[2]
        xnm = names(X)
        for(i in 1:p) {
            if(is.factor(X[[i]])) {
                Xtemp = class.ind(X[[i]])
                colnames(Xtemp) = paste(xnm[i],1:ncol(Xtemp),sep='')
                X[[i]]=Xtemp
                m=ncol(Xtemp)
                grp=c(grp, rep(m, m))
                ##grp=c(grp, rep(i, ncol(Xtemp)))
            } else {
                X[[i]]=cbind(X[[i]])
                colnames(X[[i]])=xnm[i]
                grp=c(grp, 1)
                ##grp=c(grp, i)
            }
        }
        Xtemp=cbind(X[[1]])
        if(p>1) for(i in 2:p) Xtemp=cbind(Xtemp, X[[i]])
        X=Xtemp
    }
    else if(X.class=='numeric' | X.class=='integer') {
        X=cbind(as.numeric(X))
        ##grp=1
    }
    else if(X.class=='NULL') return(X)
    else if(X.class!='matrix')
        stop('Expecting either a factor, a vector, a matrix or a data.frame')

    N <- nrow(X)
    p <- ncol(X)

    xinfo. <- matrix(nrow=p, ncol=numcut)
    nc <- numcut
    rm.vars <- c()

    if(N>0 & p>0 & (rm.const | numcut[1]>0)) {
        for(j in 1:p) {
            X.class <- class(X[1, j])[1]

            if(X.class=='numeric' | X.class=='integer') {
                xs <- unique(sort(X[ , j]))
                k <- length(xs)
                nc[j] <- numcut

                if(k %in% 0:1) {
                     rm.vars <- c(rm.vars, -j)
                     nc[j] <- 1
                     if(k==0) xs <- NA
                }
                else if(cont)
                    xs <- seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
                else if(k<numcut) {
                    xs <- 0.5*(xs[1:(k-1)]+xs[2:k])
                    nc[j] <- k-1
                }
                else if(usequants) {
                    xs <- quantile(X[ , j], type=type,
                                   probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]
                    ##xs <- quantile(X[ , j], type=type, probs=(1:numcut)/(numcut+1))
                    names(xs) <- NULL
                }
                else xs <-
                         seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
            }
            else
                stop(paste0('Variables of type ', X.class, ' are not supported'))

            ##nc[j] <- length(xs)
            xinfo.[j, 1:nc[j] ] <- xs
        }
    }

    X <- data.matrix(X)

    if(length(xinfo)>0) {
        if(is.list(xinfo)) for(j in 1:p) xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
        else if(is.matrix(xinfo)) xinfo. <- xinfo
        else stop('Only a list or a matrix can be provided for xinfo')

        for(j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
    }

    xinfo <- xinfo.

    if(rm.const && length(rm.vars)>0 && !all((1:p)==(-rm.vars))) {
        X <- X[ , rm.vars]
        nc <- nc[rm.vars]
        xinfo <- xinfo[rm.vars, ]
        grp <- grp[rm.vars]
    }
    else if(length(rm.vars)==0 || all((1:p)==(-rm.vars)))
        rm.vars <- 1:p

    dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)

    if(numcut==0) return(X)
    else return(list(X=X, numcut=as.integer(nc), rm.const=rm.vars,
                     xinfo=xinfo, grp=grp))
}
