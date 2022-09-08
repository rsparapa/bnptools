## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## predict.nft.R

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

predict.nft = function(
                       ## data
                       object,
                       x.test=object$x.train,
                       ## multi-threading
                       tc=getOption("mc.cores", 1), ##OpenMP thread count
                       ## current process fit vs. previous process fit
                       XPtr=TRUE,
                       ## predictions
                       K=0,
                       events=object$events,
                       FPD=FALSE,
                       probs=c(0.025, 0.975),
                       take.logs=TRUE,
                       na.rm=FALSE,
                       seed=NULL,
                       ## default settings for NFT:BART/HBART/DPM
                       fmu=object$fmu,
                       soffset=object$soffset,
                       drawDPM=object$drawDPM,
                       ## etc.
                       ...)
{
    attr(object, 'class') <- 'nft2'
    object$xftrain=object$x.train
    object$xstrain=object$x.train
    object$x.train=NULL
    object$xifcuts=object$xicuts
    object$xiscuts=object$xicuts
    object$xicuts=NULL
    np=length(object$x.test)
    if(np>0) {
        object$xftest=object$x.test
        object$xstest=object$x.test
        object$x.test=NULL
    }
    return(predict(
                       ## data
                       object,
                       xftest=x.test,
                       xstest=x.test,
                       ## multi-threading
                       tc=tc, ##OpenMP thread count
                       ## current process fit vs. previous process fit
                       XPtr=XPtr,
                       ## predictions
                       K=K,
                       events=events,
                       FPD=FPD,
                       probs=probs,
                       take.logs=take.logs,
                       na.rm=na.rm,
                       seed=seed,
                       ## default settings for NFT:BART/HBART/DPM
                       fmu=fmu,
                       soffset=soffset,
                       drawDPM=drawDPM,
                       ## etc.
                       ...))
}
