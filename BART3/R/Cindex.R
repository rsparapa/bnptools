## Copyright (C) 2022 Rodney A. Sparapani

## This file is part of BART3.
## Cindex.R

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

Cindex=function(risk, times, delta=NULL)
{   
    N=length(risk)
    if(N!=length(times))
        stop('risk and times must be the same length')
    if(length(delta)==0) delta=rep(1, N)
    else if(N!=length(delta))
        stop('risk and delta must be the same length')

    l=0
    k=0
    for(i in 1:N) {
        h=which((times[i]==times & delta[i]>delta) |
                (times[i]<times & delta[i]>0))
        if(length(h)>0) {
            l=l+sum(risk[i]>risk[h])
            k=k+length(h)
        }
    }
    return(l/k)
}
