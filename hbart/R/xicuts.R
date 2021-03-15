## Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
##                         Hugh A. Chipman and Rodney A. Sparapani

## This file is part of hbart.
## xicuts.R

## hbart is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.

## hbart is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author contact information
## Matthew T. Pratola: mpratola@gmail.com
## Robert E. McCulloch: robert.e.mculloch@gmail.com
## Hugh A. Chipman: hughchipman@gmail.com
## Rodney A. Sparapani: rsparapa@mcw.edu

xicuts = function(x, transposed=FALSE, numcut=100) {
    numcut.=0
    grid=list()
    if(transposed) {
        p=nrow(x)
        n=ncol(x)
        for(i in 1:p) grid[[i]]=unique(sort(x[i, ]))
        names.=dimnames(x)[[1]]
    } else {
        p=ncol(x)
        n=nrow(x)
        for(i in 1:p) grid[[i]]=unique(sort(x[ , i]))
        names.=dimnames(x)[[2]]
    }
    ##return(grid)
    if(length(names.)==0) names.=paste0('x', 1:p)
    xicuts.=list()
    for(i in 1:p) {
        numcut.[i]=length(grid[[i]])-1
        if(numcut.[i]==0)
            warning(paste0('The following column is constant:', i))
        if(numcut.[i]>=numcut || numcut.[i]==(n-1)) {
            xinc=(grid[[i]][numcut.[i]+1]-grid[[i]][1])/(numcut+1)
            xicuts.[[i]]=(1:numcut)*xinc+grid[[i]][1]
        } else {
            xicuts.[[i]]=double(numcut.[i])
            for(j in 1:numcut.[i])
                xicuts.[[i]][j]=mean(grid[[i]][c(j, j+1)])
        }
    }
    names(xicuts.)=names.
    class(xicuts.)="BARTcutinfo"
    return(xicuts.)
}

