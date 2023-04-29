## Copyright (C) 2023 Rodney A. Sparapani

## This file is part of nftbart.
## concordance.R

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

concordance=function(draws, times, delta=NULL)
{   
    N=ncol(draws)
    if(N!=length(times))
        stop('draw columns and times must be the same length')
    if(length(delta)==0) delta=rep(1, N)
    else if(N!=length(delta))
        stop('draw columns and delta must be the same length')

    C=double((N*(N-1)/2))
    ##C=matrix(nrow=N, ncol=N)
    k=1
    for(i in 1:(N-1))
        for(j in (i+1):N) {
            if((times[i]==times[j] && delta[i]>delta[j]) ||
                (times[i]<times[j] && delta[i]>0))
                C[k]=mean(draws[ , i]<draws[ , j])
            else if((times[i]==times[j] && delta[i]<delta[j]) ||
                (times[i]>times[j] && delta[j]>0))
                C[k]=mean(draws[ , i]>draws[ , j])
            else C[k]=NA
            k=k+1
        }
    return(mean(C, na.rm=TRUE))
}
