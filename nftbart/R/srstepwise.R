## Copyright (C) 2021 Rodney A. Sparapani

## This file is part of nftbart.
## nft.R

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

srstepwise=function(x, times, delta=rep(1, length(times)),
                    sle=0.15, sls=0.15, dist='lognormal',
                    type='right') {
    x=cbind(x)

    P=ncol(x)

    in.=c()
    out.=1:P
    
    step=TRUE
    H=0
    L=P

    if(dist=='gaussian') {
        dist='lognormal'
        y=Surv(time=exp(times), event=delta, type=type)
    } else {
        y=Surv(time=times, event=delta, type=type)
    }
    
    while(step) {
        fits=as.list(1:L)
        A=matrix(nrow=L, ncol=H+1)
        for(i in 1:L) {
            if(H==0)
                fits[[i]]=survreg(y~x[ , out.[i]], dist=dist)
            else 
                fits[[i]]=survreg(y~x[ , out.[i]]+x[ , in.], dist=dist)
            A[i, ]=pnorm(-abs(coef(fits[[i]])[2:(H+2)]/
                              sqrt(diag(vcov(fits[[i]]))[2:(H+2)])))
        }
##return(A)
        
        ##forward step
        j=which(order(A[ , 1])==1)[1]
        if(A[j, 1]<sle) {
            in.=c(out.[j], in.)
            H=H+1
            out.=out.[-j]
            L=L-1
            if(L==0) step=FALSE
        }
        else if(H==0) step=FALSE
        else {
            fits[[1]]=survreg(y~x[ , in.], dist=dist)
            A=A[ , -(H+1)]
            A[1, ]=pnorm(-abs(coef(fits[[i]])[2:(H+1)]/
                              sqrt(diag(vcov(fits[[i]]))[2:(H+1)])))
            j=1
            step=FALSE
        }
        
        ##backward step
        if(H>0) {
            k=which(order(-A[j, ])==1)[1]
            if(A[j, k]>=sls) {
                out.=c(out., in.[k])
                L=L+1
                in.=in.[-k]
                H=H-1
                if(step && k==1) step=FALSE
                else step=TRUE
            }
        }
    }
    
    return(in.)
}
