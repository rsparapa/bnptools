## Copyright (C) 2022 Rodney A. Sparapani

## This file is part of nftbart.
## bMM.R

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

bMM=function(X, numcut=0L, usequants=FALSE, type=7, xicuts=NULL, rm.const=FALSE,
             rm.dupe=FALSE, method="spearman", use="pairwise.complete.obs"
             ##cont=FALSE, rm.vars=NULL
             ) {
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
                grp[i]=0
                for(k in j+(1:ncol(Xtemp)))
                    ##if(!(k %in% -rm.vars))
                    grp[i]=grp[i]+1
                j=j+ncol(Xtemp)
            } else {
                X[[i]]=cbind(X[[i]])
                colnames(X[[i]])=xnm[i]
                grp[i]=1 ##-(j %in% -rm.vars)
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
    
    if(length(dimnames(X)[[2]])==0) dimnames(X)[[2]]=paste0('x', 1:ncol(X))  
    N <- nrow(X)
    p <- ncol(X)
    nc <- numcut
    ##rm.vars <- c()

    ##if(N>0 && p>0 && (rm.const || numcut[1]>0)) 
    if(N>0 && p>0 && numcut[1]>0) 
        for(j in 1:p) {
            X.class <- class(X[1, j])[1]

            if(X.class=='numeric' | X.class=='integer') {
                xs <- unique(sort(X[ , j]))
                k <- length(xs)
                nc[j] <- numcut

                if(k>1) {
                    if(k<numcut) {
                        xs <- 0.5*(xs[1:(k-1)]+xs[2:k])
                        nc[j] <- k-1
                    }
                    else if(usequants) {
                        xs <- quantile(X[ , j], type=type,
                                       probs=(0:(numcut+1))/
                                           (numcut+1))[-c(1, numcut+2)]
                        names(xs) <- NULL
                    }
                    else xs <- seq(xs[1], xs[k],
                                   length.out=numcut+2)[-c(1, numcut+2)]
                }
            }
            else
                stop(paste0('Variables of type ', X.class, ' are not supported'))

        }

    ##X <- data.matrix(X)

    if(numcut>0 && length(xicuts)==0)
        xicuts=xicuts(X, numcut=numcut)
    
    if(numcut==0) return(X)
    else {
        dummy=matrix(0, nrow=2, ncol=p)
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

        chv=NULL
        rm.vars=c()
        rm.drop=c()
        if(rm.const) {
            for(i in 1:p)
                if(nc[i]==1 & all(X[ , i]==xicuts[[i]][1])) {
                    ##warning(paste0('column ', i, ' is constant: removed'))
                    rm.vars=c(rm.vars, i)
                }
        }
        if(rm.dupe) { 
            chv=cor(X, method='spearman', use="pairwise.complete.obs")
            flag=TRUE
            while(flag) {
                count=rep(0, p)
                h=1
                same=matrix(nrow=0, ncol=2)
                for(i in 1:(p-1))
                    for(j in (i+1):p) 
                        if(!(i %in% rm.drop) && !(j %in% rm.drop) && abs(chv[i, j])==1) {
                            same=rbind(same, c(i, j))
                            count[i]=count[i]+1
                            count[j]=count[j]+1
                        }
                flag=(sum(count)>0)
                if(flag) {
                    l=which(count==max(count))[1]
                    rm.drop=c(rm.drop, l)
                }
            }
        }

        drop=unique(sort(c(rm.vars, rm.drop), TRUE))
        if(length(drop)>0) {
            X=cbind(X[ , -drop])
            nc=nc[-drop]
            chv=chv[-drop, -drop]
            dummy=dummy[ , -drop]
            for(i in drop) {
                xicuts[[i]]=NULL
                dummy[1, which(dummy[1, ]>i)]=dummy[1, which(dummy[1, ]>i)]-1
                dummy[2, which(dummy[2, ]>i)]=dummy[2, which(dummy[2, ]>i)]-1
            }
        }
        
        return(list(X=X, numcut=as.integer(nc), rm.const=sort(rm.vars),
                    xicuts=xicuts, ##grp=grp,
                    dummy=dummy, chv=chv, rm.dupe=sort(rm.drop)))
    }
}
